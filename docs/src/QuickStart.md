# Quick Start

`NMRflux.jl` brings vendor aware loading, coordinate aware data, classical NMR processing, and spin simulation into one Julia workflow. This page walks the basic path: install the package, load a dataset into `SpectData`, process a 1D FID with a composable `Chain`, and see how the same public interfaces support the external RINSE machine learning demonstrator.

## 1. Installation

```julia
import Pkg
Pkg.add(url = "https://github.com/utzgroup/NMRflux.jl.git")
```

Then:

```julia
using NMRflux
```

## 2. Load a dataset into SpectData

NMRflux reads Bruker, JEOL, Magritek Spinsolve, Varian/Agilent, and Oxford Instruments data. The example below uses the bundled Bruker dataset, so nothing needs downloading first.

```@example brukerEg
using NMRflux
using NMRflux.Examples
using Plots: plot, savefig

data_jeol = NMRflux.Examples.Data["Spheroid culture medium"]
jdf_file = joinpath(data_jeol["path"], "yp-5-fu-2.5-100.jdf")
params_JEOL, data_td = NMRflux.load(jdf_file, :JEOL)

plot(coords(data_td, 1), real.(data_td.dat);
    xlabel = "time / s",
    ylabel = "signal (a.u.)",
    title = "Bruker FID (real part)")

savefig("quickstart_bruker_fid.svg"); nothing
```

![](quickstart_bruker_fid.svg)

`NMRflux.load` returns two things: a dictionary of acquisition parameters straight from the vendor file, and a `SpectData` object holding the time domain signal together with its time axis. Which kind of path each vendor expects, and what each loader does to the raw data, is set out in [Loading NMR Data](DataLoading.md).

## 3. Apply a 1D processing chain

```@example brukerEg
N_new = max(length(data_td.dat), 2^16)

p = Chain(
    ZeroFill([N_new]),
    Apodize([0.5]),
    FourierTransform([N_new], [1]; fftshift=true),
    PhaseCorrect(-0.55pi,2pi*0.00175,1),
    MedianBaselineCorrect(1; wdw=2<<12),
)

s = p(data_td)
s ./= sqrt(sum(s.*conj(s)))
plot(
    coords(s,1) / 600 .+ 4.835,
    real(s),
    xaxis=:flip, xlabel = "1H chemical shift [ppm]", xlims=[-0.5,7.5],
    label=nothing,
    grid=false,
    yaxis=false,
    minorticks=10,
    title = "quickstart_processing_pipeline"
)

savefig("quickstart_processing_pipeline.svg"); nothing # hide
```
![](quickstart_processing_pipeline.svg)

The phase values here belong to this dataset and nothing else. For a new experiment, either read them off an interactive plot or let `AutoPhaseCorrectChen` find them. [Classical Processing Pipeline](DataProcessing.md) covers every processor in detail. The same `SpectData` representation runs through the whole framework, so the output of one processor feeds straight into the next without rebuilding axes or swapping containers.

## 4. Machine learning use case: RINSE

RINSE is a separately developed machine learning application that uses the public NMRflux interfaces from end to end. It is **not part of the NMRflux package or core API**. The RINSE workflow uses `SpinSim` for the underlying spin simulations, `SpectData` for coordinate aware signals and spectra, NMRflux processors for spectral preprocessing, and its own machine learning code for the learned restoration model. That separation is deliberate. NMRflux supplies the light weight reusable framework; RINSE shows one way to build a specialised application on top of it. See [RINSE](RINSE.md) for the synthetic data generation workflow, the dataset format, the training and inference procedure, and the Zenodo dataset citation.
