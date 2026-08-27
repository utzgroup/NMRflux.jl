###########################################################################
# File:        DataSet.jl
# Project:     NMRflux.jl
#
# Description:
#   Defines the SpectData array type and associated coordinate handling.
#   Provides the high level load interface for importing NMR data from
#   supported vendor formats into SpectData objects.
#
# Original implementation:
#   Marcel Utz
#
# Additional vendor support:
#   Manaz Kaleel
###########################################################################

# abstract type CoordMap{T} <: AbstractVector{T} end

@doc raw"""
    struct SpectData{T,N} <: AbstractArray{T,N}

Basic  structure for spectral data.
"""
struct SpectData{T,N} <: AbstractArray{T,N}
    dat::AbstractArray{T,N}
    coord::NTuple{N,AbstractVector}
end

import Base.size
import Base.getindex
import Base.setindex!
import Base.IndexStyle
import Base.showarg

size(S::SpectData) = size(S.dat)
getindex(S::SpectData, k::Integer) = getindex(S.dat,k)
setindex!(S::SpectData, v, k::Integer) = setindex!(S.dat,v,k)
IndexStyle(S::SpectData) = IndexStyle(S.dat)

Base.showarg(io::IO, A::SpectData, toplevel) = print(io, typeof(A), " with coords:", A.coord)

## support broadcasting ------------------------------------------------------
Base.BroadcastStyle(::Type{<:SpectData}) = Broadcast.ArrayStyle{SpectData}()

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SpectData}}, ::Type{ElType}) where ElType
    # Scan the inputs for the SpectData:
    A = find_spdta(bc)
    # Use the size and coordinates of the SpectData to construct a new one with the appropriate element type:
    SpectData(similar(Array{ElType}, axes(bc)), A.coord)
end

"`A = find_spdta(As)` returns the first SpectData among the arguments."
find_spdta(bc::Base.Broadcast.Broadcasted) = find_spdta(bc.args)
find_spdta(args::Tuple) = find_spdta(find_spdta(args[1]), Base.tail(args))
find_spdta(x) = x
find_spdta(::Tuple{}) = nothing
find_spdta(a::SpectData, rest) = a
find_spdta(::Any, rest) = find_spdta(rest)

## ----------------------------------------------------------------------------

## support slicing ------------------------------------------------------------
getindex(S::SpectData{T,N},I::Vararg{<:Integer,N}) where {T,N} = getindex(S.dat,CartesianIndex(I...))

function getindex(S::SpectData,I...)
    adat=S.dat[I...]
    if all(x->x isa Integer,I) || all(x->x isa CartesianIndex,I)
        return adat
    else
        newcoords=[]
        for (k,i) in enumerate(I)
            if !(i isa Integer) # integer indices are dropped
                push!(newcoords,S.coord[k][i])
            end
        end
    
    return SpectData(adat,(newcoords...,))
    end
end 

## ----------------------------------------------------------------------------

@doc raw"""
    function coords(S::SpectData)
    function coords(S::SpectData,k::Integer)

returns a tuple with the coordinates of `S`, analogous to `axes(S)`.
The second form returns the coordinate of the `k`-th dimension.
"""
coords(S::SpectData) = S.coord
coords(S::SpectData,k::Integer) = S.coord[k]

function SpectData(A::AbstractArray)
    sz=size(A)
    coord=map(x->1:x,sz)
    return SpectData(A,coord)
end

import Base: convert,similar
convert(::Type{SpectData},A::AbstractArray) = SpectData(A)


@doc raw"""
    function load(path::String,vendor::Symbol)

Load a data set located at `path`, and return a dictionary with the 
acquisition parameters as well as a `SpectData` object with the actual
(time domain) data. `vendor` designates
the origin data format. Currently implemented are

- `:Bruker`: the path points to a directory with a Bruker NMR data set.
- `:JEOL`: the path points to a JEOL `.jdf` file
- `:Magritek`: the path points to a Magritek `data.1d` and `acqu.par` files
- `:Varian`: the path points to a Varian `fid` and `procpar` files
- `:Oxford`: the path points to a Oxford `.dx/jdx` files
"""
function load(f::String,vendor::Symbol)
    if vendor == :Bruker 
        params = FileIO.readBrukerParameterFile(f*"/acqus")
        rawdata = FileIO.readBrukerFID(f*"/fid")
        rawdata = rawdata[params["GRPDLY"]:end]
        tcoord = range(0.0,step=1.0/params["SW_h"],length=length(rawdata))
        return params,SpectData(rawdata,(tcoord,))
    elseif vendor == :JEOL
        io = open(f,"r")
        header,params,data = FileIO.readJEOL(io)
        close(io)
        arr, coords = FileIO.reshapeJEOL(header, params, data)
        return params, SpectData(arr, coords)
    elseif vendor == :Magritek
        parameterFile = joinpath(f, "acqu.par")
        dataFile = joinpath(f, "data.1d")

        params = FileIO.readMagritekParameterFile(parameterFile)
        haskey(params, "bandwidth") || error("Magritek acqu.par is missing the " * "bandwidth parameter")

        _header, _storedAxis, rawdata = FileIO.readMagritekFID(dataFile)
        spectralWidth = Float64(params["bandwidth"]) * FileIO.KILOHERTZ_TO_HERTZ

        spectralWidth > 0 || error("Magritek bandwidth must be positive")
        tcoord = range(0.0, step=1.0 / spectralWidth, length=length(rawdata),)

        return params, SpectData(rawdata, (tcoord,))
    elseif vendor == :Varian
        parameterFile = joinpath(f, "procpar")
        dataFile = joinpath(f, "fid")

        params = FileIO.readVarianParameterFile(parameterFile)

        haskey(params, "sw") ||
            error(
                "Varian procpar is missing the spectral width sw"
            )

        _fileHeader, _blockHeader, rawdata =
            FileIO.readVarianFID(dataFile)

        spectralWidth = Float64(params["sw"])

        spectralWidth > 0 ||
            error("Varian spectral width sw must be positive")

        if haskey(params, "np")
            procparPointCount =
                Int(round(Float64(params["np"])))

            binaryPointCount =
                2 * length(rawdata)

            if procparPointCount != binaryPointCount
                @warn(
                    "Varian procpar np differs from the binary file",
                    procparPointCount,
                    binaryPointCount,
                )
            end
        end

        tcoord = range(
            0.0,
            step=1.0 / spectralWidth,
            length=length(rawdata),
        )

        return params, SpectData(rawdata, (tcoord,))
    elseif vendor == :Oxford
        dataFile =
            FileIO.resolveOxfordJcampFile(f)

        records, timeAxis, rawdata =
            FileIO.readOxfordFID(dataFile)

        params =
            FileIO.oxfordParameterDictionary(records)

        return params, SpectData(rawdata, (timeAxis,))
    else 
        error("Unsupported data format")
    end
end

