# Roadmap to NMRflux.jl 1.0

`NMRflux.jl` aims to be a stable, vendor neutral framework for NMR data representation, processing, and spin dynamics simulation, with interfaces that stay easy to combine with external scientific and machine learning tools.

## Completed core capabilities

- **Coordinate aware data**
  - `SpectData` subtypes `AbstractArray`.
  - Numerical values and physical coordinates stay associated.
  - Standard Julia indexing, slicing, and broadcasting all work.

- **Vendor neutral data access**
  - High level loading for Bruker, JEOL, Magritek Spinsolve, Varian/Agilent, and Oxford Instruments data.
  - Vendor readers stay reachable through `NMRflux.FileIO`.
  - Supported datasets convert into the common `SpectData` representation.

- **Composable processing**
  - `NMRProcessor` and `NMRProcessor1D` provide the processor abstractions.
  - `Chain` assembles processing operations into reusable pipelines.
  - Core operations cover zero filling, apodization, Fourier transformation, manual and automatic phase correction, baseline correction, differentiation, cumulative integration, and peak alignment.

- **Spin dynamics simulation**
  - `SpinSim` provides spin operators, coupling operators, density operator propagation, FID simulation, transition resolved spectra, and Lorentzian line shapes.
  - `SpinSim.FID` and `SpinSim.PeakSpect` return `SpectData`, so simulated data feeds the same processing pipeline as experimental data. `SpinSim.RungeKutta` and `SpinSim.Propagate` return plain arrays of sampled observables, which the caller pairs with an acquisition time axis. `SpinSim.Spectrum` returns transition frequencies and amplitudes, which `PeakSpect` renders onto an axis.
  - `GISSMO` builds `SpinSim` Hamiltonians from parametrised spin systems in the GISSMO database.

- **Julia ecosystem interoperability**
  - `SpectData` follows the array interface instead of introducing a separate machine learning data model.
  - NMRflux outputs pass to external Julia numerical and machine learning libraries through ordinary array operations.

## Applied examples and illustrative workflows

Application specific workflows can build on the public NMRflux interfaces without becoming part of the core package. **RINSE restoration demonstrator:** RINSE is a separately developed machine learning application that uses the NMRflux simulation, data representation, and processing interfaces. Its application specific scripts and its `GenerateFIDs` utility ship with RINSE and define nothing in the NMRflux API. The workflow is documented in [RINSE](RINSE.md).

## In progress

- **Documentation and examples**
  - Docstrings for the exported processors that still lack them, and for the exported `SpinSim` routines.
  - Curated example datasets for reproducible tutorials.
  - Doctests across the public API.

- **Public API review**
  - Review the exported names and naming conventions.
  - Check consistency between the top level processing interfaces and the `SpinSim`, `FileIO`, `GISSMO`, and `Examples` modules.
  - Settle the sign and phase conventions shared by the propagation and transition routes in `SpinSim`.
  - Verify the two dimensional JEOL section ordering against a real dataset.
  - Stabilise the interfaces before 1.0.

## Planned release timeline

| Version | Planned features | Target release |
|---:|---|---|
| **0.1.0** | Core cleanup, basic documentation, processing examples | Q3 2026 |
| **0.2.0** | Refined `SpectData`, transforms, plotting integration | Q4 2026 |
| **0.3.0** | Extended vendor and metadata support, loading helpers | Q1 2027 |
| **0.4.0** | Interoperability examples, including an external ML demonstrator | Q2 2027 |
| **1.0.0** | Stable public interfaces, complete manual and API reference | TBD |
