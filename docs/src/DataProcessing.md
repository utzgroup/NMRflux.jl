# Classical Processing Pipeline

Processing operations in `NMRflux.jl` are callable objects, not plain functions. Each carries its own configuration and acts on `SpectData`. Once configured, they compose into pipelines you can reuse across datasets.

## 1. The NMRProcessor abstraction

```julia
abstract type NMRProcessor <: Function end
```

A concrete processor stores whatever parameters or precomputed plans it needs, and defines a call method that uses them. It behaves like a function, but carries its configuration with it. Processors that act on a single dimension subtype `NMRProcessor1D` instead. These only have to handle a one dimensional `SpectData`, and they must carry a `dim` field naming the dimension they apply to. `mapslices` does the rest, lifting them to higher dimensional data. A processor will also accept a bare `AbstractArray`: the array is promoted to `SpectData` with unit range coordinates, processed, and the numerical part handed back.

## 2. Chaining processors

`Chain` applies its arguments in the order you list them:

```julia
Chain(fs::Vararg{Function}) = reduce(∘, reverse(fs))
```

So `Chain(p1, p2, p3)` behaves like `x |> p1 |> p2 |> p3`. The signature accepts any `Function`, so plain functions and anonymous functions can sit in a chain alongside processors. We use the bundled JEOL example throughout this page.

```@example JEOLEg
using NMRflux
using NMRflux.Examples
using Plots: plot, savefig

data_jeol = NMRflux.Examples.Data["Spheroid culture medium"]
jdf_file = joinpath(data_jeol["path"], "yp-5-fu-2.5-100.jdf")
params_JEOL, data_td_JEOL = NMRflux.load(jdf_file, :JEOL)

size(data_td_JEOL.dat)
```

```@example JEOLEg
plot(coords(data_td_JEOL, 1), real.(data_td_JEOL.dat);
    xlabel = "time / s",
    ylabel = "signal (a.u.)",
    title = "JEOL FID (real part)"
    )

savefig("loaded_FID.svg"); nothing # hide
```
![](loaded_FID.svg)

## 3. Zero filling

`ZeroFill([SI...])` extends the data by appending zeros, one target size per dimension. It raises the digital resolution of the transformed spectrum without adding any experimental information. The target size must be at least the current size.

```@example JEOLEg
N_orig = length(data_td_JEOL.dat)
N_new = max(N_orig, 2^16)

zf = ZeroFill([N_new])
data_td_JEOL_zf = zf(data_td_JEOL)

(size(data_td_JEOL.dat), size(data_td_JEOL_zf.dat))
```

Passing `:` for a dimension keeps its current length. The processor resolves that `:` to a concrete size the first time you call it and remembers the result, so build a fresh `ZeroFill` for each dataset of a different length.

## 4. Apodization

`Apodize` multiplies the time domain signal by an exponential weighting function,

```math
f(t) = \exp(-Rt),
```

where `R` is the decay rate along the selected dimension. If the coordinate is in seconds, `R` is in inverse seconds. Pass one entry per dimension, and use `:` to leave a dimension alone: `Apodize([0.5, :])` weights only the first dimension of a two dimensional dataset.

!!! warning "In place modification"
    `Apodize` writes the weighted signal back into the array it was handed, and the object it returns shares that same array. The input is therefore modified in place. Keep a copy of the raw FID if you need it later, and do not apply the same processor twice to the same data.

```@example JEOLEg
ap = Apodize([0.5])
data_td_JEOL_zf_ap = ap(data_td_JEOL_zf)

plot(coords(data_td_JEOL_zf_ap, 1), real.(data_td_JEOL_zf_ap.dat);
    xlabel = "time / s",
    ylabel = "signal (a.u.)",
    title = "JEOL FID after (ZF + AP)"
    )

savefig("JEOL_td_zf_ap_plot.svg"); nothing # hide
```
![](JEOL_td_zf_ap_plot.svg)

## 5. Fourier transform

`FourierTransform(SI, dims; fftshift=true)` transforms along the listed dimensions and replaces each transformed coordinate with a frequency coordinate. The sizes in `SI` are used to build an FFTW plan once, which makes repeated transforms of same sized datasets cheap. Build a new processor when the size changes.

The new coordinate runs from `-SW/2` to `+SW/2`, where `SW` is the reciprocal of the original sample spacing. Zero frequency therefore sits in the centre of the axis, which is what `fftshift=true` arranges in the data. Keep the default: the coordinate is centred either way, so turning the shift off leaves the axis labelling out of step with the data.

```@example JEOLEg
SI = [length(data_td_JEOL_zf_ap.dat)]
ft = FourierTransform(SI, [1]; fftshift=true)

data_fd_JEOL_zf_ap = ft(data_td_JEOL_zf_ap)
size(data_fd_JEOL_zf_ap.dat)
```

```@example JEOLEg
plot(coords(data_fd_JEOL_zf_ap, 1), real.(data_fd_JEOL_zf_ap.dat);
    xaxis = :flip,
    xlabel = "frequency / Hz",
    ylabel = "signal (a.u.)",
    title = "JEOL spectrum (ZF + AP + FT)"
    )

savefig("JEOL_fd_zf_ap_plot.svg"); nothing # hide
```
![](JEOL_fd_zf_ap_plot.svg)

### 5.1. PPM conversion

The frequency axis that `FourierTransform` produces is centred on the transmitter. Dividing it by the observe frequency in MHz converts it to ppm, and adding the shift of the reference line puts it on the conventional scale. This dataset was acquired at 600 MHz with the transmitter on the water resonance at 4.835 ppm, which is where the two constants below come from.

```@example JEOLEg
s = data_fd_JEOL_zf_ap
s ./= sqrt(sum(s.*conj(s)))

plot(
    coords(s,1) / 600 .+ 4.835,
    real(s),
    xaxis=:flip, xlabel = "1H chemical shift [ppm]", xlims=[-0.5,7.5],
    label=nothing,
    grid=false,
    yaxis=false,
    minorticks=10,
    title = "JEOL spectrum (ZF + AP + FT)"
    )

savefig("JEOL_fd_zf_ap_plot_ppm.svg"); nothing # hide
```
![](JEOL_fd_zf_ap_plot_ppm.svg)

## 6. Phase correction

`PhaseCorrect(ph0, ph1, dim)` multiplies the spectrum by `exp(i*ph0) * exp(i*ph1*f)`, where `f` is the frequency coordinate of dimension `dim`. `ph0` is the zero order phase in radians. `ph1` is the first order coefficient, so its units are radians per unit of the coordinate, which for a spectrum in Hz means radians per Hz.

```@example JEOLEg
ph0 = -0.55pi
ph1 = 2pi * 0.00175

pc = PhaseCorrect(ph0, ph1, 1)
data_fd_JEOL_zf_ap_pc = pc(data_fd_JEOL_zf_ap)

s = data_fd_JEOL_zf_ap_pc
s ./= sqrt(sum(s.*conj(s)))

plot(
    coords(s,1) / 600 .+ 4.835,
    real(s),
    xaxis=:flip, xlabel = "1H chemical shift [ppm]", xlims=[-0.5,7.5],
    label=nothing,
    grid=false,
    yaxis=false,
    minorticks=10,
    title = "JEOL spectrum (ZF + AP + FT + PC)"
    )

savefig("JEOL_fd_zf_ap_pc_plot.svg"); nothing # hide
```
![](JEOL_fd_zf_ap_pc_plot.svg)

## 7. Baseline correction

`MedianBaselineCorrect(dim; wdw=4096, stp=32)` estimates a slowly varying baseline from the local extrema of the real part, smooths it with a Gaussian kernel, and subtracts it. The method follows Friedrichs, *J. Biomol. NMR* **5** (1995), 147-153. `wdw` is the half width in points of the median window, and also the half width at which the Gaussian kernel is truncated. The standard deviation of that kernel is `wdw/sqrt(50)`, roughly a seventh of the truncation width. `stp` is accepted and stored but the current implementation does not use it.

The processor keeps the real part and drops the imaginary one, so its output is real valued. Run it after any step that still needs the complex spectrum.

```@example JEOLEg
mbc = MedianBaselineCorrect(1; wdw=2048)
data_fd_JEOL_zf_ap_pc_bc = mbc(data_fd_JEOL_zf_ap_pc)

size(data_fd_JEOL_zf_ap_pc_bc.dat)
```

```@example JEOLEg
s = data_fd_JEOL_zf_ap_pc_bc
s ./= sqrt(sum(s.*conj(s)))

plot(
    coords(s,1) / 600 .+ 4.835,
    real(s),
    xaxis=:flip, xlabel = "1H chemical shift [ppm]", xlims=[-0.5,7.5],
    label=nothing,
    grid=false,
    yaxis=false,
    minorticks=10,
    title = "JEOL spectrum (ZF + AP + FT + PC + BC)"
    )

savefig("JEOL_fd_zf_ap_pc_bc_plot.svg"); nothing # hide
```
![](JEOL_fd_zf_ap_pc_bc_plot.svg)

## 8. Automatic phase correction

`AutoPhaseCorrectChen(dim; verbose=false, γ=1.0e-5)` estimates both phase parameters by minimising the entropy of the first derivative of the real part, following Chen et al., *J. Magn. Reson.* **158** (2002), 164-168. It runs a one dimensional search over `ph0` first, then refines both parameters together. The `γ` parameter adds a penalty on the negative values of the phase corrected derivative, which helps in noisy spectra. Setting `γ=0.0` switches the penalty off. The method works best when it only has to refine an approximate manual correction, so apply a rough `PhaseCorrect` first.

```@example JEOLEgAPC
using NMRflux
using NMRflux.Examples
using Plots: plot, savefig

data_jeol = NMRflux.Examples.Data["Spheroid culture medium"]
jdf_file = joinpath(data_jeol["path"], "yp-5-fu-2.5-100.jdf")
params_JEOL, data_td_JEOL = NMRflux.load(jdf_file, :JEOL)

# Zero filling
N_orig = length(data_td_JEOL.dat)
N_new = max(N_orig, 2^16)
zf = ZeroFill([N_new])
data_td_JEOL_zf = zf(data_td_JEOL)

# Apodization
ap = Apodize([0.5])
data_td_JEOL_zf_ap = ap(data_td_JEOL_zf)

# FFT
SI = [length(data_td_JEOL_zf_ap.dat)]
ft = FourierTransform(SI, [1]; fftshift=true)
data_fd_JEOL_zf_ap = ft(data_td_JEOL_zf_ap)

# Auto phase correction
pc_guess = PhaseCorrect(-0.55pi, 2pi * 0.00175, 1)
data_pc_guess = pc_guess(data_fd_JEOL_zf_ap)

# the manual guess above is already close and this spectrum has good SNR,
# so the entropy term alone is enough here
apc = AutoPhaseCorrectChen(1; verbose=false, γ=0.0)

data_apc = apc(data_pc_guess)

# PPM conversion
s = data_apc
s ./= sqrt(sum(s.*conj(s)))

plot(
    coords(s,1) / 600 .+ 4.835,
    real(s),
    xaxis=:flip, xlabel = "1H chemical shift [ppm]", xlims=[-0.5,7.5],
    label=nothing,
    grid=false,
    yaxis=false,
    minorticks=10,
    title = "JEOL spectrum (ZF + AP + FT + PC + AutoPC)"
    )

savefig("JEOL_fd_zf_ap_apc_plot.svg"); nothing # hide
```
![](JEOL_fd_zf_ap_apc_plot.svg)

Always look at the result. Strong baseline distortion, a large residual solvent signal, or unusual line shapes can all pull the optimiser to a solution that scores well on entropy and looks wrong on screen.

## 9. Other 1D processors

Three further processors are exported and follow the same calling pattern. All three run below against the baseline corrected spectrum from section 7, so the documentation build exercises them every time it runs.

- `Derivative(dim)` returns the first derivative along `dim` from a five point central difference stencil. The first two and last two points fall outside the stencil, so they come back as edge artefacts.

- `Integral(dim)` returns the running integral along `dim`, accumulated point by point with the rectangle rule. The result is a cumulative trace the same length as the input. Getting a peak area from it means taking the difference between two points, so treat it as raw material for an integration workflow.

Both read the sample spacing with `step`, so the coordinate of that dimension has to be a range.

```@example JEOLEg
deriv = Derivative(1)(data_fd_JEOL_zf_ap_pc_bc)
cumulative = Integral(1)(data_fd_JEOL_zf_ap_pc_bc)

(size(deriv.dat), size(cumulative.dat))
```

- `PeakAlign(dim, readpos, wdw)` searches the `wdw` points either side of `readpos`, takes the largest absolute value it finds there, and shifts the spectrum cyclically so that feature lands on `readpos`. This is the usual way to bring a set of spectra onto a common reference signal such as TMS.

The search window is not clipped at the array bounds, so keep `readpos` at least `wdw` points away from either end of the spectrum. The shift is computed in whole points, so the aligned feature lands within one point of `readpos`.

```@example JEOLEg
aligned = PeakAlign(1, 0.0, 500)(data_fd_JEOL_zf_ap_pc_bc)
size(aligned.dat)
```

## 10. A typical 1D pipeline

The common sequence is zero filling, apodization, Fourier transformation, phase correction, and baseline correction. Assembled as a chain, that reads:

```@example JEOLEgPipeline
using NMRflux
using NMRflux.Examples
using Plots: plot, savefig

data_jeol = NMRflux.Examples.Data["Spheroid culture medium"]
jdf_file = joinpath(data_jeol["path"], "yp-5-fu-2.5-100.jdf")
params_JEOL, data_td_JEOL = NMRflux.load(jdf_file, :JEOL)

N_new = max(length(data_td_JEOL.dat), 2^16)

p = Chain(
    ZeroFill([N_new]),
    Apodize([0.5]),
    FourierTransform([N_new], [1]; fftshift=true),
    PhaseCorrect(-0.55pi, 2pi * 0.00175, 1),
    MedianBaselineCorrect(1; wdw=2048),
)

JEOL_data_processed = p(data_td_JEOL)

s = JEOL_data_processed
s ./= sqrt(sum(s.*conj(s)))

plot(
    coords(s,1) / 600 .+ 4.835,
    real(s),
    xaxis=:flip, xlabel = "1H chemical shift [ppm]", xlims=[-0.5,7.5],
    label=nothing,
    grid=false,
    yaxis=false,
    minorticks=10,
    title = "JEOL spectrum after the processing chain"
    )

savefig("JEOL_full_pipeline_plot.svg"); nothing # hide
```
![](JEOL_full_pipeline_plot.svg)

The chain reproduces the step by step result from the sections above. Once a set of parameters works for a given experiment, the whole pipeline becomes a single reusable object.
