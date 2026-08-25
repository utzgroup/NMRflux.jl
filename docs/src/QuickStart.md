# QuickStart
`NMRflux.jl` provides a unified Julia workflow for NMR data loading, classical processing, spin based simulation, and machine learning (Flux). This QuickStart shows the "happy path":

1. Install `NMRflux.jl`
2. Load an example dataset into `SpectData`
3. Apply a typical 1D processing pipeline with `Chain`
4. Run a toy deep learning denoiser using `SpectData` + Flux

# 1. Installation
To install the package from GitHub:

```@julia
import Pkg
Pkg.add(url = "https://github.com/utzgroup/NMRflux.jl.git")
```

Once installation is complete:

```@julia
using NMRflux
```

# 2. Load a dataset into SpectData

`SpectData` is the central data structure in `NMRflux.jl`. It stores:
- `sd.dat` : the numerical array (time domain FID or spectrum)
- `sd.coord`: coordinate vectors for each dimension (time, frequency, ppm, etc)

The documentation ships with small example datasets accessible via `NMRflux.Examples.Data`.

```@example brukerEg
using NMRflux
using NMRflux.Examples
using Plots: plot, plot!, savefig

data_bruker = NMRflux.Examples.Data["HCC cell culture media spectra"]

# High level vendor loader (recommended)
params_bruker, data_td = NMRflux.load(joinpath(data_bruker["path"], "10"), :Bruker)

t = data_td.coord[1]
y = real.(data_td.dat)

plot(t, y;
xlabel = "time / s",
ylabel = "signal (a.u.)",
title = "Bruker FID (real part)")

savefig("quickstart_bruker_fid.svg"); nothing # hide
```
![](quickstart_bruker_fid.svg)

# 3. Classical 1D processing pipeline

In `NMRflux.jl`, processing is implemented via `NMRProcessor` functors that can be composed using `Chain`.

Here is a minimal example of a processing pipeline:

```@example brukerEg

Processing = Chain(
    ZeroFill([2^16]),
    FourierTransform([2^16],[1]),
    AutoPhaseCorrectChen(1)
)

proc = Processing(data_td / 1e10 )

plot(proc.coord[1]/700 .+ 4.78, real.(proc),xaxis=:flip,
xlims=[-0.5,9.0]
)

savefig("quickstart_basic_spect.svg") ; nothing # hide
```
![](quickstart_basic_spect.svg)

`Processing` consists of zero filling to 64k points, Fourier transformation,
and automatic phase correction.

In practice, other steps may be added, for example a digital filter to 
remove the solvent artefact (appearing in the middle of the spectrum),
apodization (to balance resolution and sensitivity), and baseline 
correction (to remove distortions due to probe ringing):

```@example brukerEg
dt = step(data_td.coord[1])

Processing = Chain(
    ZeroFill([2^16]),
    Apodize([0.5π]),
    DigitalFilter(NMRflux.BandReject(-0.0025,0.005,1024),1),
    FourierTransform([2^16],[1]),
    PhaseCorrect(0.0,2pi*1024*dt,1),
    AutoPhaseCorrectChen(1),
    MedianBaselineCorrect(1,wdw=2048)
)

proc = Processing(data_td / 1e10 )

integral = proc |> Integral(1)

plot(proc.coord[1]/700 .+ 4.78, real.(proc),xaxis=:flip, xlims=[-0.5,9.0],label="Spectrum", xlabel="Chemical Shift (ppm)")
plot!(proc.coord[1]/700 .+ 4.78, real(integral ) ./20 , xaxis=:flip,label="Integral")

savefig("quickstart_processing_pipeline.svg"); nothing # hide
```
![](quickstart_processing_pipeline.svg)

# 4. Toy deep learning example (SpectData + Flux)
We generate a target spectrum using NMRflux's classical baseline correction, then train a small Flux model to approximate this mapping. This demonstrates seamless interoperability between `SpectData`, `NMRflux` processing pipelines, and `Flux` models. This toy example demonstrates end-to-end compatibility between:
- `SpectData` (data container)
- `NMRflux` processing tools (FFT pipeline)
- `Flux` (a minimal 1D conv denoiser)

```@julia
using NMRflux, NMRflux.Examples, Flux, Statistics


d=NMRflux.Examples.Data["HCC cell culture media spectra"] # Shipped Bruker example dataset
_, td = NMRflux.load(joinpath(d["path"], "10"), :Bruker)  # Load time domain FID as SpectData

N = max(length(td.dat), 2^16)                            # Zero fill to at least 64k points                                      
sd = NMRflux.Chain(                                       # Standard 1D NMR processing pipeline                                                       
    NMRflux.ZeroFill([N]),
    NMRflux.Apodize([0.5]),
    NMRflux.FourierTransform([N], [1]; fftshift=true)
)(td)

y = Float32.(real.(sd.dat))                                            # Raw data
yt = Float32.(real.(NMRflux.MedianBaselineCorrect(1; wdw=256)(sd).dat)) # Target generation

x = Float32.(collect(eachindex(y)) ./ length(y))  # Normalised frequency coordinate
m = Flux.Chain(Dense(1,16,tanh), Dense(16,1))     # Minimal MLP mapping position -> signal
opt = Adam(1e-2)                                  # Adam optimiser with small learning rate

for _ in 1:200
    loss() = Flux.mse(vec(m(reshape(x,1,:))), yt) # MSE between prediction and target
    gs = Flux.gradient(loss, Flux.params(m))      # Compute gradients of the loss w.r.t. model
    Flux.update!(opt, Flux.params(m), gs)         # Update model parameters using Adam
end
```
# End of Manual
