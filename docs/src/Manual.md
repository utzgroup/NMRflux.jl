# User Manual

`NMRflux.jl` provides a common framework for vendor aware NMR data loading, coordinate aware data representation, composable processing, and spin dynamics simulation.This manual maps the framework and points at the page that covers each piece in full.

## 1. Installation

```julia
import Pkg
Pkg.add(url = "https://github.com/utzgroup/NMRflux.jl.git")

using NMRflux
```

## 2. Data loading

`NMRflux.load(path, vendor)` converts a supported vendor dataset into `SpectData` and hands back the acquisition parameters alongside it. The current loading paths cover Bruker, JEOL, Magritek Spinsolve, Varian/Agilent, and Oxford Instruments data. What `path` has to point to differs by vendor, so check the table on the loading page before your first attempt.

The `NMRflux.FileIO` submodule holds the vendor readers themselves and stays available whenever the high level path hides something you need. See [Loading NMR Data](DataLoading.md).

## 3. SpectData

`SpectData` is the coordinate aware data representation shared by everything else. It behaves like a Julia array while carrying the physical coordinate for each dimension, and `coords` is the exported accessor for those coordinates.

See [Working with SpectData](SpectData.md).

## 4. Processing

Processing operations subtype `NMRProcessor`, or `NMRProcessor1D` when they act on a single dimension. They are callable objects that compose with `Chain`. The package ships zero filling, apodization, Fourier transformation, manual and automatic phase correction, baseline correction, differentiation, cumulative integration, and peak alignment.

See [Classical Processing Pipeline](DataProcessing.md).

## 5. Spin dynamics

`SpinSim` carries the spin dynamics functionality: spin operators, coupling operators, time domain propagation, FID simulation, and transition resolved spectra. The `GISSMO` module sits next to it and turns entries from the GISSMO database into `SpinSim` Hamiltonians.

See [Spin Dynamics](SpinDynamics.md).

## 6. Example data

The `NMRflux.Examples` module ships the datasets used throughout this documentation. `NMRflux.Examples.Data` is a dictionary keyed by dataset name, and each entry records the dataset directory, the files it contains, and any tags from the accompanying TOML description. Run `keys(NMRflux.Examples.Data)` to see what is available in your installation.

## 7. Machine learning interoperability

NMRflux does not ship a machine learning model, and it does not need one. Its array oriented interfaces are built to interoperate with external Julia tooling, so NMR data can be handed to a learning library without a second data model in between. RINSE is included in this documentation as a worked application level demonstrator. RINSE and its synthetic data generator are distributed separately from NMRflux and form no part of the core API.

See [RINSE](RINSE.md).
