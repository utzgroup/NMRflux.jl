module NMRflux

    export SpectData,coords,load
    export NMRProcessor, Chain, FourierTransform, Apodize, ZeroFill,
           PhaseCorrect, MedianBaselineCorrect, Derivative, Integral,
           AutoPhaseCorrectChen, PeakAlign, DigitalFilter, CoordMap
   
    export SpinSim
    export GISSMO
    export Craft

    include("DataSet.jl")
    include("Examples.jl")
    include("NMRProcessor.jl")
    include("FileIO.jl")
    include("SpinSim.jl")
    include("GISSMO.jl")
    include("Craft.jl")
end
