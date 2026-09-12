module FluxExt

using NMRflux
using Flux

function Flux.gpu(A::SpectData)
    return SpectData(Flux.gpu(A.dat), A.coord)
end

function Flux.cpu(A::SpectData)
    return SpectData(Flux.cpu(A.dat), A.coord)
end

end # module    