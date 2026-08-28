module PlotsExt

using NMRflux
using Plots

const SpectData1D = NMRflux.SpectData{<:Number,1}

Plots.plot(A::SpectData1D, args...; kwargs...) =
    Plots.plot(NMRflux.coords(A, 1), real.(A.dat), args...; kwargs...)

Plots.plot!(A::SpectData1D, args...; kwargs...) =
    Plots.plot!(NMRflux.coords(A, 1), real.(A.dat), args...; kwargs...)

Plots.plot!(plt::Plots.Plot, A::SpectData1D, args...; kwargs...) =
    Plots.plot!(plt, NMRflux.coords(A, 1), real.(A.dat), args...; kwargs...)

end
