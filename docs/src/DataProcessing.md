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

N_new = max(N_orig, N_target)       # Never shrink: only zero fill if N_orig < N_target
zf = ZeroFill([N_new])              # Create ZeroFill processor

data_p = data_td_joel |> zf  # Apply zero filling to the SpectData object
(size(data_td_joel), size(data_p)) # Before/after sizes
```

```@example joelEg
t_zf = data_p.coord[1]
y_zf = real.(data_p.dat)

plot(t_zf, y_zf;
xlabel = "time / s",
ylabel = "signal (a.u.)",
title = "JOEL FID after zero filling")

savefig("joel_fid_zf_plot.svg"); nothing  # hide
```
![](joel_fid_zf_plot.svg)

## 1.4 Apodization (`Apodize`)
Apodization applies a decay function (window) to the time-domain FID. In practice, this damps the tail of the FID, reducing truncation artefacts and high frequency noise in the frequency domain at the cost of some line broadening. It is typically applied *after* zero filling and before the Fourier transform. 

In `NMRflux.jl`, apodization is implemented by the processor Apodize, which
acts on `SpectData` by multiplies the data along selected dimensions by an exponential factor of the form:
```@julia
f(t) = exp.(-R * t)
```
where:
- The `t` is the coordinate axis of the processed dimension (usually time)
- The `R` is the user specified decay constant for that dimension

Internally, Apodize uses the coordinate vector of each selected dimension to compute this exponential weighting and multiplies it into the underlying data array. The coordinate vectors are preserved. For time domain FIDs, `R` therefore has units of `1/seconds`. Continuing from the previous section, we apply apodization to the zero filled FID `data_td_joel_zf`:
```@example joelEg
ap = Apodize([0.5]) # Decay constant for the first (time) dimension

data_p = data_td_joel |> zf |> ap ;

```
This produces a windowed time domain signal suitable for Fourier transformation.

```@example joelEg
# Extract time axis and real part AFTER apodization
t_ap = data_p.coord[1]
y_ap = real.(data_p)

plot(t_ap, y_ap; xlabel="time / s", ylabel="signal (a.u.)", title="JOEL FID after ZF + AP")
savefig("joel_fid_zf_ap_plot.svg"); nothing  # hide
```
![](joel_fid_zf_ap_plot.svg)

## 1.5 Fourier transform (FourierTransform)
The Fourier transform converts a time domain FID into a frequency domain spectrum. After zero filling and apodization, applying the FFT produces a complex spectrum whose real and imaginary parts can be used for further processing (phase correction, baseline correction, peak picking, etc.).

In `NMRflux.jl`, the processor FourierTransform wraps FFTW's FFT planning and
updates the coordinate axes accordingly. For each transformed dimension:
- The FFT is applied to the data
- The coordinate axis is replaced by a frequency axis based on the sampling interval (Nyquist theorem)
- By default, the spectrum is shifted so that zero frequency appears at the centre (fftshift = true). This behaviour can be disabled by setting fftshift = false.

The `FourierTransform` constructor is declared as:

```@julia
function FourierTransform(SI::Vector, dims::Vector; fftshift = true)
     dummy = zeros(ComplexF64, SI...)
     plan = FFTW.plan_fft(dummy, dims)
     return FourierTransform(dims, SI, fftshift, plan)
end
```
- `SI`: size of the data array (e.g. [N] for 1D, [N1, N2] for 2D)
- `dims`: dimensions along which the FFT is computed (e.g. [1], [1, 2])
- `fftshift`: whether to apply FFTW.fftshift so that zero frequency is in the centre of the axis

Example: Continuing from Section *4.4*, we start from the apodized, zero filled FID `data_td_joel_zf_ap`:
```@example joelEg
SI = 2^16           # Size of the apodized time domain data (1D)
ft = FourierTransform([SI], [1]; fftshift = true) # Construct a FourierTransform along the first dimension, with fftshift

data_p = data_td_joel  |> zf |> ap |> ft  ;   # Apply FT to the apodized zero filled SpectData
```

```@example joelEg
f_ap = data_p.coord[1]   # frequency axis (Hz)
y_ap = real.(data_p) # real spectrum

plot(f_ap, y_ap, xaxis=:flip,
xlabel = "frequency [Hz]",
ylabel = "signal (a.u.)",
title = "JOEL spectrum (ZF + AP + FT)")

savefig("joel_fd_zf_ap_plot.svg"); nothing # hide
```
![](joel_fd_zf_ap_plot.svg)

## 1.6 Phase correction (PhaseCorrect)
After Fourier transformation, NMR spectra usually require phase correction before the real part of the spectrum can be interpreted as an absorption mode spectrum. Phase errors can lead to dispersive lineshapes, asymmetric peaks, and negative peak components. `PhaseCorrect` applies zero-order and first-order phase correction along a chosen dimension.

- `ph0`: zero order phase (radians), a uniform rotation applied to all points
- `ph1`: first order phase (radians per axis unit), a linear phase ramp
- `dim`: the dimension along which the correction is applied (typically `1` for 1D spectra)

Continuing from the previous section, we start from the frequency domain, zero filled, apodized spectrum `data_fd_joel_zf_ap`:
```@example joelEg
ph0 = -0.55pi       # zero-order phase (radians)
ph1 = 2pi*0.00175   # first-order phase (radians)
dim = 1             # apply along the first (frequency) dimension

pc = PhaseCorrect(ph0, ph1, dim)
data_p =data_td_joel  |> zf |> ap |> ft |> pc 

f_pc = data_p.coord[1]
y_pc = real.(data_p)

plot(f_pc, y_pc, xaxis=:flip,
     xlabel = "frequency [Hz]",
     ylabel = "signal (a.u.)",
     title  = "JOEL spectrum (ZF + AP + FT + PC)")

savefig("joel_fd_zf_ap_pc_plot.svg"); nothing # hide
```
![](joel_fd_zf_ap_pc_plot.svg)


In practice, `ph0` and `ph1` are adjusted until the peaks in the real part of the spectrum are predominantly absorptive and symmetric. This can be done manually, interactively, or by using an automatic phase correction method. Via the generic `NMRProcessor1D` machinery, `PhaseCorrect` can also be applied slice wise along a chosen dimension of higher dimensional SpectData objects.

## 1.7 Automatic phase correction (`AutoPhaseCorrectChen`)
AutoPhaseCorrectChen provides an automatic phase correction processor based on the minimum entropy method of Chen et al. [J. Magn. Reson. 158 (2002), 164–168]. The method estimates phase correction parameters by minimizing an entropy based objective function.

The processor first optimizes the zero-order phase and then refines both zero-order and first-order phase terms. The optional penalty parameter `γ` adds a penalty for negative real valued points in the corrected spectrum. This can help reduce overcorrection, especially in noisy spectra.

- `dim`: dimension along which phase correction is applied, typically `1` for 1D spectra
- `verbose`: if `true`, prints optimisation information
- `γ`: penalty strength for negative real-valued points; default is `1.0e-5`

Automatic phase correction is most reliable when the spectrum is already reasonably close to the correct phase. For spectra with substantial first-order phase error, an approximate manual correction can be applied first, followed by `AutoPhaseCorrectChen` as a refinement step.

Example: automatic refinement after approximate manual phase correction.

We start again from the zero filled, apodized, Fourier transformed spectrum `data_fd_joel_zf_ap`. First, an approximate manual correction is applied. The automatic method is then used to refine the result:

```@example joelEg
ph0_guess = -0.55pi
ph1_guess = 2pi*0.00175

pc_guess = PhaseCorrect(ph0_guess, ph1_guess, 1)

apc = NMRflux.AutoPhaseCorrectChen(1, verbose=false, γ=0.0e-5)
data_p = data_td_joel |> zf |> ap |> ft |> pc_guess |> apc ;

f_apc = data_p.coord[1]
y_apc = real.(data_p)

plot(f_apc, y_apc, xaxis=:flip,
     xlabel = "frequency [Hz]",
     ylabel = "signal (a.u.)",
     title  = "JOEL spectrum (ZF + AP + FT + PC + automatic PC refinement)")

savefig("joel_fd_zf_ap_apc_plot.svg"); nothing # hide
```
![](joel_fd_zf_ap_apc_plot.svg)

The automatically corrected spectrum should still be inspected visually. This is especially important for noisy spectra, spectra with strong baseline distortions, spectra with large residual solvent signals, or spectra with unusual peak shapes. In such cases, changing the penalty parameter `γ` may help reduce excessive negative peaks or overcorrection.

For example, a stronger penalty can be used as follows:

```@example joelEg
apc_stronger_penalty = NMRflux.AutoPhaseCorrectChen(1, verbose=false, γ=1.0e-4)
data_fd_joel_zf_ap_apc_penalty = apc_stronger_penalty(data_p)

eltype(data_fd_joel_zf_ap_apc_penalty.dat), ndims(data_fd_joel_zf_ap_apc_penalty.dat)
```

## 1.7 Baseline correction (MedianBaselineCorrect)
After phase correction, spectra often exhibit slowly varying offsets or slopes in the real part, known as baseline distortions. These can bias peak integration and make peak picking less reliable. Baseline correction should be applied after phase correction, so that the real part contains the absorptive peaks. The `MedianBaselineCorrect` implements a robust baseline correction for the real part of a spectrum. The algorithm follows the method of M. S. Friedrichs (Journal of Biomolecular NMR, 5 (1995) 147-153) and proceeds as follows along a chosen dimension:
- For each position along the axis, a local window of width `2 * wdw + 1` points is considered
- Within this window, local extrema (minima and maxima) of the real part are extracted
- The median of these extrema is taken as a robust estimate of the local baseline
- The sequence of baseline estimates is then smoothed by convolution with a Gaussian like kernel
- The resulting smooth baseline is subtracted from the original complex data (affecting primarily the real part)

The processor is constructed as:
```@julia
NMRflux.MedianBaselineCorrect(dim; wdw=2<<12)
```
where:
- `dim`: dimension along which baseline correction is performed (typically 1 for 1D spectra)
- `wdw`: half width of the local window in points (controls the "smoothness" scale of the baseline)

Example: Continuing from the previous section, we start from the phased frequency domain spectrum `data_fd_joel_zf_ap_pc`:
```@example joelEg
mbc = NMRflux.MedianBaselineCorrect(1; wdw=2<<12) # baseline correction along frequency dimension
data_p = data_td_joel |> zf |> ap |> ft |> pc |> apc |> mbc;

```

```@example joelEg
f_bc = data_p.coord[1]          # frequency axis (Hz)
y_bc = real.(data_p)        # real part after baseline correction

plot(f_bc, y_bc,  xaxis=:flip,
     xlabel = "frequency [Hz]",
     ylabel = "signal (a.u.)",
     title  = "JOEL spectrum (ZF + AP + FT + PC + BC)")

savefig("joel_fd_zf_ap_pc_bc_plot.svg"); nothing # hide
```
![](joel_fd_zf_ap_pc_bc_plot.svg)

Here:
- `dim = 1` selects the first dimension (the frequency axis) for baseline correction
- `wdw = 256` sets the half width of the local window; larger values produce a smoother, more slowly varying baseline estimate

The output `data_fd_jeol_zf_ap_pc_bc` is a `SpectData` object with the same coordinates as the input but with the baseline of the real part largely removed, providing a cleaner spectrum for peak picking and integration.

## 1.8 A typical 1D processing pipeline
In practice, the processing steps described above are often combined into a single pipeline that transforms a raw time domain FID into a clean, baseline corrected spectrum ready for peak picking and integration.

For 1D data, a typical sequence is:
1. Zero filling
2. Apodization
3. Fourier transform
4. Phase correction
4. Baseline correction

The `NMRflux.jl` allows these processors to be combined using `Chain`, which applies them in sequence.

```@example joelEg
# (1) Loading JOEL data
data_joel = NMRflux.Examples.Data["Spheroid culture medium"]["files"][1]
params_joel, data_td_joel = NMRflux.load(data_joel,:JEOL) ;

# (1) Zero filling
N_orig   = length(data_td_joel.dat)              # Original number of points
N_target = 2^16                                  # Target size: 64k points
N_new    = max(N_orig, N_target)                 # Never shrink: only zero fill if N_orig < N_target
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
