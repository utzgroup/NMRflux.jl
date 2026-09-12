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
NMRflux.FourierTransformPlan
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

## Craft API
```@autodocs
Modules = [NMRflux.Craft]
Order   = [:module, :type, :function, :method]
```

## Index

```@index
```
