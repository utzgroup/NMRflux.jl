# Reference

```@meta
CurrentModule = NMRflux
```

## Core NMRflux API

```@docs
NMRflux.SpectData
NMRflux.coords
NMRflux.load
NMRflux.NMRProcessor
NMRflux.NMRProcessor1D
NMRflux.Chain
NMRflux.FourierTransform
NMRflux.MedianBaselineCorrect
NMRflux.AutoPhaseCorrectChen
NMRflux.Derivative
NMRflux.Integral
NMRflux.PeakAlign
NMRflux.conv
NMRflux.entropy
NMRflux.extrema
NMRflux.find_spdta
```

### Processors without docstrings

Three exported processors carry no docstring in the source yet, so they cannot appear in the block above. Their signatures are:

```julia
ZeroFill(SI::Vector{Union{Integer,Colon}})
Apodize(R::Vector{Union{Real,Colon}})
PhaseCorrect(ph0::Float64, ph1::Float64, dim::Int32)
```

`ZeroFill` takes one target size per dimension and pads with zeros; `:` keeps a dimension at its current length. `Apodize` takes one exponential decay rate per dimension and multiplies the data by `exp(-R*t)` along each of them; `:` leaves a dimension untouched. `PhaseCorrect` multiplies dimension `dim` by `exp(i*ph0) * exp(i*ph1*f)`, where `f` is that dimension's coordinate. [Classical Processing Pipeline](DataProcessing.md) works through all three.

## FileIO API

The vendor readers are exposed through `NMRflux.FileIO`. The parsing helpers that support them are deliberately left out of this reference.

```@docs
NMRflux.FileIO
NMRflux.FileIO.readBrukerFID
NMRflux.FileIO.readBrukerParameterFile
NMRflux.FileIO.readJEOL
NMRflux.FileIO.reshapeJEOL
NMRflux.FileIO.readMagritekFID
NMRflux.FileIO.readMagritekParameterFile
NMRflux.FileIO.readVarianFID
NMRflux.FileIO.readVarianParameterFile
NMRflux.FileIO.readOxfordBlocks
NMRflux.FileIO.readOxfordFID
NMRflux.FileIO.resolveOxfordJcampFile
NMRflux.FileIO.oxfordParameterDictionary
```

## SpinSim API

```@docs
NMRflux.SpinSim
NMRflux.SpinSim.Kron
NMRflux.SpinSim.Spectrum
NMRflux.SpinSim.PeakSpect
NMRflux.SpinSim.clorentzian
NMRflux.SpinSim.expm
NMRflux.SpinSim.FID
```

### SpinSim routines without docstrings

The remaining exported names in `SpinSim` carry no docstring in the source yet. `SpinOp`, `TwoSpinOp`, `OpJstrong`, and `OpJweak` build operators, with `⊗` as the infix alias for `Kron`; `Commutator` and `Trc` are the linear algebra helpers behind the propagators; `RungeKutta` and `Propagate` are the two time domain routines; and `Sx`, `Sy`, `Sz`, `Sp`, `Sm`, and `Id` are the single spin matrices. [Spin Dynamics](SpinDynamics.md) covers all of these apart from `Commutator` and `Trc`.

## GISSMO API

```@docs
NMRflux.GISSMO
NMRflux.GISSMO.Hamiltonian
NMRflux.GISSMO.SpinMatrix
NMRflux.GISSMO.search
NMRflux.GISSMO.NMRsignals
```

## Examples API

```@docs
NMRflux.Examples
NMRflux.Examples.Data
```

## Index

```@index
```
