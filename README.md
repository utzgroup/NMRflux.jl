# NMRflux.jl

[![Build Status](https://github.com/utzgroup/NMRflux.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/utzgroup/NMRflux.jl/actions/workflows/CI.yml?query=branch%3Amain)

**NMRflux.jl** is a Julia framework for processing, simulating, and machine-learning-assisted
analysis of nuclear magnetic resonance (NMR) data. It provides a unified, vendor-neutral
interface that covers the full workflow from raw instrument files to publication-ready spectra
and ML-ready datasets.

---

## Motivation

NMR spectroscopy generates large, complex datasets that require a chain of non-trivial
processing steps before any scientific information can be extracted. In quantitative or
high-throughput settings — such as metabolomics, metabolic flux analysis, or microfluidic NMR —
this chain must be reproducible, composable, and straightforward to extend with custom
processing logic or machine learning models.

Existing NMR software is largely vendor-specific, GUI-driven, and not designed for
programmatic or automated workflows. Julia's scientific computing ecosystem, by contrast,
offers high performance, expressive type systems, and seamless interoperability with
modern ML libraries. NMRflux.jl was built to bring these advantages to NMR data science.

The package grew from concrete research needs in the Utz group (KIT), consolidating tools
that had been developed piecemeal for spin simulation, spectral processing, and deep
learning-based denoising into a single, consistent framework. It supersedes an earlier
internal toolkit (NMRlab.jl / NMR.jl) and is designed to be useful both inside and outside
the group.

---

## Design Principles

### 1. A single array-like data container: `SpectData`

All data — raw FIDs, zero-filled arrays, frequency-domain spectra, 2D datasets — is
represented by the `SpectData{T,N}` type:

```julia
struct SpectData{T,N} <: AbstractArray{T,N}
    dat   :: AbstractArray{T,N}    # the numerical data (complex or real)
    coord :: NTuple{N,AbstractVector}  # one coordinate vector per dimension
end
```

Because `SpectData` subtypes `AbstractArray`, it works transparently with broadcasting,
slicing, plotting libraries, and Flux without any adapter code. Physical axes (time in
seconds, frequency in Hz, ppm) travel with the data through every processing step,
eliminating the bookkeeping errors that arise when arrays and their axes are managed
separately.

### 2. Composable processing via `NMRProcessor` functors

Every processing operation is a callable struct that subtypes `NMRProcessor <: Function`.
Functors carry their parameters and any pre-computed state (e.g. FFTW plans) as struct
fields. Individual processors are assembled into pipelines with `Chain`:

```julia
pipeline = Chain(ZeroFill([N]), Apodize([0.5]), FourierTransform([N], [1]), PhaseCorrect(ph0, ph1, 1))
spectrum = pipeline(fid)
```

`Chain` applies processors left-to-right, equivalent to `fid |> p1 |> p2 |> ...`. This
design keeps pipelines readable, testable step-by-step, and easy to extend with new
processor types.

### 3. Vendor-neutral data loading

A single high-level entry point, `NMRflux.load(path, :Bruker)` or
`NMRflux.load(path, :JEOL)`, reads instrument-specific binary files, applies necessary
corrections (e.g. Bruker group-delay compensation), constructs the correct physical time
axis, and returns a ready-to-use `SpectData` object alongside a parameter dictionary.
Low-level readers in `NMRflux.FileIO` are also available for users who need direct access
to raw arrays and acquisition parameters.

### 4. Spin dynamics simulation in Hilbert space

The `SpinSim` submodule implements NMR spin evolution in Hilbert space using the density
operator formalism. Hamiltonians are built from sparse matrix representations of spin
operators, making the computation tractable for systems of moderate size while retaining
physical accuracy. `SpinSim` underpins both the GISSMO metabolite interface and the
synthetic FID generator.

### 5. Integrated ML workflows

`SpectData` objects can be passed directly to Flux models. A `GenerateFIDs` submodule
produces paired clean/corrupted FIDs from a TOML configuration file, providing an
effectively unlimited source of supervised training data for spectral denoising,
artefact correction, and related ML tasks. The package also ships with curated example
datasets (Bruker and JEOL) that are available immediately after installation.

---

## Getting Started

### Installation

NMRflux.jl requires Julia 1.6 or later. Install from GitHub:

```julia
import Pkg
Pkg.add(url = "https://github.com/utzgroup/NMRflux.jl.git")
```

Then load the package:

```julia
using NMRflux
```

### Load and inspect a Bruker FID

```julia
using NMRflux, NMRflux.Examples

# Bundled example dataset
d = NMRflux.Examples.Data["HCC cell culture media spectra"]

# High-level loader: returns acquisition parameters and a SpectData FID
params, fid = NMRflux.load(joinpath(d["path"], "10"), :Bruker)

t = fid.coord[1]        # time axis in seconds
y = real.(fid.dat)      # real part of the FID
```

### Process a 1D spectrum

```julia
using Plots

N = max(length(fid.dat), 2^16)

pipeline = Chain(
    ZeroFill([N]),
    Apodize([0.5]),
    FourierTransform([N], [1]; fftshift = true),
    PhaseCorrect(0.0, 0.0, 1),
    NMRflux.MedianBaselineCorrect(1; wdw = 256)
)

spectrum = pipeline(fid)

plot(spectrum.coord[1], real.(spectrum.dat),
     xaxis = :flip, xlabel = "frequency / Hz", ylabel = "signal (a.u.)")
```

Each processor in the chain updates both the data array and the coordinate axes, so
`spectrum.coord[1]` automatically holds the correctly scaled frequency axis after the
Fourier transform.

For a JEOL `.jdf` file the workflow is identical — only the load call changes:

```julia
params, fid = NMRflux.load("path/to/spectrum.jdf", :JEOL)
```

Full documentation, including tutorials for spin simulation, GISSMO queries, and the
ML denoising pipeline, is available in `docs/src/`.

---

## Features

### Currently implemented

**Data container**
- `SpectData{T,N} <: AbstractArray{T,N}` — N-dimensional data type carrying physical
  coordinate axes alongside the numerical array.
- Native support for broadcasting, slicing, and integration with Julia ecosystem tools
  (Plots, Flux, IJulia, Pluto).

**Vendor-neutral file I/O**
- High-level loader (`NMRflux.load`) for Bruker TopSpin (FID + acqus) and JEOL (.jdf)
  formats, returning `SpectData` with a correctly scaled time axis.
- Low-level readers in `NMRflux.FileIO` for direct access to raw binary data and
  acquisition parameter dictionaries (supporting both 32- and 64-bit Bruker formats).

**Classical processing pipeline**
- `ZeroFill` — pads the FID to a target length and extends the time axis consistently.
- `Apodize` — applies exponential line broadening along selected dimensions.
- `FourierTransform` — FFTW-backed FFT with optional fftshift; updates frequency axes
  automatically.
- `PhaseCorrect` — zero- and first-order phase correction.
- `AutoPhaseCorrectChen` — automatic phase correction by entropy minimisation
  (Chen et al., *J. Magn. Reson.* 2002).
- `MedianBaselineCorrect` — robust local-median baseline subtraction
  (Friedrichs, *J. Biomol. NMR* 1995).
- `Derivative`, `Integral`, `PeakAlign` — additional spectral manipulation tools.
- `Chain` — composes any sequence of processors into a single callable pipeline.

**Spin dynamics simulation (`SpinSim`)**
- Hilbert-space density operator simulation using sparse matrix Hamiltonians.
- Construction of spin Hamiltonians from chemical shifts and J-coupling constants.
- Frequency-domain and time-domain spectrum computation.

**GISSMO metabolite interface (`GISSMO`)**
- Queries the GISSMO spin-system database (BMRB) by compound ID over HTTP.
- Parses XML spin-simulation files and constructs the corresponding sparse Hamiltonians
  for direct use with `SpinSim`.

**Synthetic FID generation (`GenerateFIDs`)**
- TOML-configured batch generation of synthetic 1H NMR FIDs.
- Produces paired clean / artefact-corrupted FIDs for supervised ML training.
- Artefacts include phase errors, solvent residual signals, baseline distortions, and
  controlled Gaussian noise.
- Output is a 2D `SpectData` object, directly compatible with downstream processing
  and Flux-based models.

**Machine learning integration**
- `SpectData` is natively compatible with Flux models — no adapter layer required.
- Helpers for batching and array-format conversion (WHCN, Cartesian, FFT-domain).
- End-to-end examples combining `NMRflux` processing pipelines with Flux training loops.

**Example datasets**
- Bundled Bruker (HCC cell culture media) and JEOL (spheroid culture medium) datasets
  accessible via `NMRflux.Examples.Data`.

---

### Planned and in progress

The package is currently at version 1.0.0-DEV. The following features are either
partially implemented or planned for upcoming releases:

| Version | Focus area |
| ------: | ---------- |
| **0.2** | Refined `SpectData` types; improved axis handling; extended plotting helpers |
| **0.3** | Additional vendor formats; downloadable curated example datasets; richer metadata support |
| **0.4** | ML processing modules (spectral denoiser, FID cleaner) with documented training scripts |
| **1.0** | Stable public API; complete user manual and API reference; registered in the Julia General registry |

Specific items under active development:

- Downloadable curated FID and spectrum datasets for reproducible tutorials.
- Final API review and naming-convention unification across `NMRProcessor`, `SpinSim`,
  `FileIO`, and ML submodules.
- Extended doctest and unit-test coverage for interface stability.
- `SpecCleaner` module for automated spectral artefact removal.
- Automatic phase and baseline correction integrated into a one-call "clean spectrum"
  convenience function.

---

## Contributing and Contact

NMRflux.jl is under active development. Bug reports, feature requests, and pull requests
are welcome via the GitHub repository. For direct enquiries, contact
`marcel.utz@kit.edu`.

If you use NMRflux.jl in published work, please acknowledge this; a citable reference
will be added here once available.
