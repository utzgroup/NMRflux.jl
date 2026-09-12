# NMRflux.jl

*A modular Julia framework for NMR data handling, processing, and simulation*

`NMRflux.jl` provides a common Julia framework for handling, processing, simulating, and interpreting NMR data. The package grew out of practical need for a julia package for NMR spectroscopy: reliable vendor aware data loading, coordinate aware data representation, composable processing pipelines, and spin dynamics simulation. The same array oriented interfaces also make it straightforward to hand NMR data to external Julia libraries, machine learning packages included, without introducing a second data model along the way.

## Package features

- Vendor neutral loading of Bruker, JEOL, Magritek Spinsolve, Varian/Agilent, and Oxford Instruments data into the common `SpectData` representation.
- A coordinate aware `SpectData` type that subtypes `AbstractArray` and keeps numerical data together with its physical axes.
- Composable processing tools covering zero filling, apodization, Fourier transformation, manual and automatic phase correction, baseline correction, differentiation, cumulative integration, and peak alignment..
- Spin dynamics simulation through the general purpose `SpinSim` module, with `GISSMO` on top of it for parametrised spin systems from the GISSMO database.
- Extensible processing abstractions designed for custom algorithms and for integration with external machine learning workflows.
- A worked machine learning use case, RINSE, developed and distributed separately, showing how the public NMRflux interfaces support a complete learning workflow.

!!! note "Scope and applications"
    NMRflux.jl provides the general data, processing, and simulation framework. Application specific workflows can build on these public interfaces without becoming part of the core package. RINSE appears in this documentation as one such demonstrator. It is developed and distributed separately from NMRflux.jl.


## Manual outline

- Start with **Getting Started** for installation, vendor data loading, and a basic 1D processing pipeline.
- **User Manual** gives a short guided overview of the package, and the pages under **Advanced topics** cover data loading, `SpectData`, processing, and spin dynamics in more detail.
- **Development** holds the roadmap to 1.0, which tracks what is finished and what is still moving.
- **Reference** lists the documented public functions, types, and modules.
- **Machine learning demonstrator** documents RINSE as a use case of the public NMRflux interfaces.

## Feedback

`NMRflux.jl` is under active development. Send bug reports, feature requests, and suggestions through the project GitHub repository, or contact `marcel.utz@kit.edu`.

## Citing NMRflux.jl

If you use `NMRflux.jl` in published work, please cite the corresponding NMRflux publication once the formal reference is available.
