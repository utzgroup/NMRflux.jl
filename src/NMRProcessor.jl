
@doc raw"""
    abstract type NMRProcessor <: Function end

Abstract data type for NMR processing routines. 
New processing tools should be declared as subtypes of
`NMRProcessor`.
"""
abstract type NMRProcessor <: Function end

@doc raw"""
    (m::NMRProcessor)(A::AbstractArray)

here is some fallback behaviour. NMR processors are defined
for SpectData objects. However, they can also be applied to
any AbstractArray, by promoting it to a SpectData, and then
extracting the data part.
"""
(m::NMRProcessor)(A::AbstractArray) = m(SpectData(A)).dat


@doc raw"""
    Chain(fs::Vararg{NMRProcessor}) 

returns a chain of processing tools, which will be
applied in order (the first in the argument list is 
applied first)
"""
Chain(fs::Vararg{Function}) = reduce(∘, reverse(fs))


# TODO: Clean up API according to the following principles:
# 1. Argument structure of NMRProcessors:
#     - Dimension should be a keyword argument `dim` with default value 1
#     - Other parameters should be keyword arguments with default values, if possible
# 2. Whenever possible, functionality should be provided as NMRProcessor1D (as opposed to the more general NMRProcessor).

import FFTW

struct FourierTransform <: NMRProcessor
    dims::Vector{Integer}
    SI::Vector{Integer}
    fftshift::Bool
    plan
end

@doc raw"""
    function FourierTransform(SI::Vector,dims::Vector; fftshift=true)

Fourier transform processor for data sets of size `SI`. `dims` is a Vector
of the dimensions along which a Fourier transform will be computed.
The corresponding coordinates are automatically replaced by frequencies,
based on the Nyqvist theorem. The zero frequency appears in the centre of the
spectrum. 

The function produces a FFTW plan for the Fourier transform, which is stored ^in
the returned `FourierTransform` object. This makes it more efficient to apply
the Fourier transform to multiple data sets of the same size. If the size of the
data set changes, a new `FourierTransform` object should be created.

"""
function FourierTransform(SI::Vector,dims::Vector; fftshift=true)
    dummy=zeros(ComplexF64,SI...)
    plan=FFTW.plan_fft(dummy,dims)
    return( FourierTransform(dims,SI,fftshift,plan))
end

function (ft::FourierTransform)(S::SpectData)
    ftdat = ft.plan*S.dat
    newcoord=[]
    for (k,d) in enumerate(S.coord)
        if k in ft.dims && d isa AbstractRange
            Δf = 1.0/step(d)
            push!(newcoord, range(-Δf/2,Δf/2,length=length(d)))
        else
            push!(newcoord,d)
        end
    end
    if ft.fftshift
        ftdat=FFTW.fftshift(ftdat,ft.dims)
    end
    return SpectData(ftdat,(newcoord...,))
end


struct ZeroFill <: NMRProcessor
    SI::Vector{Union{Integer,Colon}}
end

function (zf::ZeroFill)(A::SpectData{T,N}) where {T,N}
    oldsize = size(A)
    newsize = zf.SI
    newcoord=[]
    for (k,d) in enumerate(newsize)
        if d isa Colon
            newsize[k]=oldsize[k]
            push!(newcoord,A.coord[k])
        elseif A.coord[k] isa AbstractRange
            extrange = range(first(A.coord[k]),step=step(A.coord[k]),length=d)
            push!(newcoord,extrange)
        else
            push!(newcoord,A.coord[k])
        end
    end
    newA = zeros(T, newsize...)
    oldRange = map(x->1:x,oldsize)
    newA[oldRange...] = A
    return SpectData(newA,(newcoord...,)) 
end

struct Apodize <: NMRProcessor
    R::Vector{Union{Real,Colon}}
end

function (ap::Apodize)(A::SpectData)
    apo = A.dat
    for (k,r) in enumerate(ap.R)
        if !(r isa Colon)
            ix=ones(Int64,ndims(apo))
            ix[k] = length(A.coord[k])
            f = exp.(-ap.R[k] .* A.coord[k])
            apo .*= reshape(f,ix...)
        end
    end

    return(SpectData(apo,A.coord))
end


@doc raw"""
    function CoordMap(f::Function; dim::Integer=1)

returns a processor that replaces the `dim`-th coordinate of a `SpectData` by
`f.(coord)`, i.e., applies `f` to each element of that coordinate vector. The
underlying data is left unchanged. This is useful, e.g., to rescale or relabel
an axis, such as converting a frequency axis from Hz to ppm.

**Example:**
```julia
hz_to_ppm = CoordMap(f -> f/600.13, 1)   # convert a Hz axis to ppm at 600.13 MHz
spectrum_ppm = hz_to_ppm(spectrum)
```
"""
struct CoordMap <: NMRProcessor
    f::Function
    dim::Int64
end

function CoordMap(f::Function; dim::Integer=1)
    return CoordMap(f,dim)
end

function (cm::CoordMap)(A::SpectData{T,N}) where {T,N}
    newcoord = ntuple(k -> k == cm.dim ? cm.f.(A.coord[k]) : A.coord[k], N)
    return SpectData(A.dat, newcoord)
end


@doc raw"""
    abstract type NMRProcessor1D <: NMRProcessor

Data type for processing tools that apply to a single dimension, i.e., that are inherently
1D. They need to be defined as functors that act on a `AbstractVector{T}`. They must
contain a field `.dim` that indicates which dimension in a multidimensional array they
should be applied to.
"""
abstract type NMRProcessor1D <: NMRProcessor end

function (np1d::NMRProcessor1D)(A::SpectData{T,N}) where {T,N}
    return mapslices(np1d, A; dims=np1d.dim)
end

@doc raw"""
    function FFT(; dim::Integer=1)

returns a processor that computes the Fourier transform of a spectrum along
the dimension `dim` (default `1`), via `FFTW.fft` followed by `FFTW.fftshift`
so that zero frequency appears in the centre. The coordinate along `dim` is
replaced by a frequency axis, based on the Nyqvist theorem: for a coordinate
with `n` points and step `dt`, the new coordinate runs from `-1/(2dt)` to
`1/(2dt)`.

Unlike `FourierTransform`, `FFT` requires no pre-built plan and no target
size — as an `NMRProcessor1D`, it applies to a `SpectData` of any
dimensionality via `mapslices`, transforming only `dim` and leaving all other
dimensions (data and coordinates) untouched.
"""
struct FFT <: NMRProcessor1D
    dim::Int64
end

function FFT(; dim::Integer=1)
    return FFT(dim)
end

function (fft::FFT)(A::SpectData{T,1}) where {T<:Number}
    Δf = 1.0/step(A.coord[1])
    newcoord = range(-Δf/2,Δf/2,length=length(A.coord[1]))
    return SpectData(FFTW.fftshift(FFTW.fft(A.dat)),(newcoord,))
end


@doc raw"""
    function PhaseCorrect(; ph0::Real=0.0, ph1::Real=0.0, dim::Integer=1)

returns a processor that applies zero- and first-order phase correction to a
spectrum along the dimension `dim` (default `1`), multiplying it by
`exp(im*ph0) * exp(im*ph1*f)`, where `f` is the coordinate of that dimension.
With the default `ph0=ph1=0.0`, the processor is the identity.

- `ph0`: zero order phase, in radians — a uniform rotation applied to every point.
- `ph1`: first order phase, in radians per unit of the coordinate (e.g. radians
  per Hz for a frequency-domain spectrum) — a linear phase ramp across the axis.
"""
struct PhaseCorrect <: NMRProcessor1D
    ph0::Float64
    ph1::Float64
    dim::Int32
end

function PhaseCorrect(; ph0::Real=0.0, ph1::Real=0.0, dim::Integer=1)
    return PhaseCorrect(Float64(ph0), Float64(ph1), dim)
end

function (pc::PhaseCorrect)(x::SpectData{T,1}) where {T<:Number}
    c = exp(im*pc.ph0).* exp.(im*pc.ph1.*coords(x,1))
    return x.*c
end

@doc raw"""
    function MedianBaselineCorrect(; dim::Integer=1, wdw::Integer=4096, stp::Integer=32)

returns a processor that subtracts a slowly varying baseline from the real
part of a spectrum along the dimension `dim` (default `1`), following the
algorithm of M. S. Friedrichs, *Journal of Biomolecular NMR*, **5** (1995)
147-153.

- `wdw`: half width, in points, of the local window used to estimate the
  baseline from local extrema; also the half width at which the Gaussian
  smoothing kernel is truncated.
- `stp`: accepted and stored, but not currently used by the implementation.
"""
struct MedianBaselineCorrect <: NMRProcessor1D
    dim::Int64
    wdw::Int64
    stp::Int64
    gauss::Vector{Float64}
end

function MedianBaselineCorrect(; dim::Integer=1, wdw::Integer=4096, stp::Integer=32)
    g=exp.(-25*((-wdw:wdw)./wdw).^2)
    g=g/sum(g)
    return MedianBaselineCorrect(dim,wdw,stp,g)
end

@doc """
    function extrema(X::AbstractArray{T,N}, dim::Integer) 

returns an array of booleans indicating all extremal values of `X` along the dimension `dim`.
"""
function extrema(X::AbstractArray{T,N}, dim::Integer) where {T,N} 
    shifter=zeros(Int64,ndims(X))
    shifter[dim]=1
    left=circshift(X,-shifter)
    right=circshift(X,shifter)

    minmaxima = (X .> left .&& X .> right) .|| (X .< left .&& X .< right)
    return minmaxima
end

wrap(n,l)=[mod(k,l) for k in n]


@doc"""
    function conv(X::AbstractArray{T1,N}, y::AbstractVector{T2},dim::Integer) where {N>1,T1,T2}

computes the convolution of the array `X` with the vector `y` along the dimension `dim`.
The ends of `X` are zero-padded such that the result is guaranteed to have the same size as `X`.
`NMRflux.conv()` uses a direct algorithm for the convolution, not fft. It is therefore efficient
when the length of `y` is much less than the corresponding dimension of `X`. If this is not
the case and performance is critical, a different algorithm should be used.
"""
function conv(X::AbstractArray{T1,N}, y::AbstractVector{T2},dim::Integer) where {N,T1,T2}
    return mapslices(a->conv(a,y), X, dims=dim)
end

function conv(x::AbstractVector{T1}, y::AbstractVector{T2} ) where {T1,T2}
    N=length(x)
    M=length(y)>>1
    lrange = isodd(length(y)) ? (-M:M) : (-M:(M-1))
    b = [ sum( ((k+l>0 && k+l<=N) ? x[k+l] : last(x) ) * y[l+M+1] for l=lrange)   for k=1:N]
end

import Statistics

@doc raw"""
    function (mb::MedianBaselineCorrect)(s::SpectData)

subtract baseline for the real part of `s` by the algorithm of M. S. Friedrichs,
*Journal of Biomolecular NMR*,  **5** (1995) 147  153.
"""
function (mb::MedianBaselineCorrect)(s::SpectData{T,1}) where {T<:Number}
    r=real.(s.dat)
    xtr = extrema(r,1)
    xind = findall(xtr)  # find the indices of all extrema in the spectrum
    bl=[Statistics.median(r[filter(x-> (x>=k-mb.wdw) && (x<=k+mb.wdw), xind)]) for k in 1:length(s.dat)] # find the index of the extremum closest to the 32500th point (the artifact)
    c=conv(bl,mb.gauss,1)
    return SpectData(r.-c, s.coord)
end


@doc raw"""
    function Derivative(; dim::Integer=1)

returns a processor that computes the first derivative of a spectrum along the dimension `dim` (default `1`).
"""
struct Derivative <: NMRProcessor1D
    dim::Int64
    Derivative(; dim::Integer=1) = new(dim)
end

function (der::Derivative)(spect::SpectData{T,1}) where {T<:Number}
    s=spect.dat
    x=spect.coord[1]
    dx = diff(x)
    inc = [dx[1]; dx]   # reuse the first spacing at the leading edge, instead of wrapping around
    d = 1.0/12*(8*[s[2:end];0]-8*[0;s[1:(end-1)]] - [s[3:end];0;0] + [0;0;s[1:(end-2)]] ) ./ inc
    return SpectData(d, spect.coord)
end


@doc raw"""
    function Integral(; dim::Integer=1)

returns a processor that computes the integral of a spectrum along the dimension `dim` (default `1`).
"""
struct Integral <: NMRProcessor1D
    dim::Int64
    Integral(; dim::Integer=1) = new(dim)
end

function (int::Integral)(spect::SpectData{T,1}) where {T<:Number}
    s=spect.dat
    x=spect.coord[1]
    dx = diff(x)
    inc = [dx[1]; dx]   # reuse the first spacing at the leading edge, instead of wrapping around
    d = cumsum(s .* inc)
    return SpectData(d, spect.coord)
end


ent(x) = -x*log(x)

import Optim

@doc raw"""
    function `entropy(s::SpectData{T,1})` 
        
computes the entropy of the first derivative
in the real part of an
NMR spectrum as defined by Chen et al. in
*Journal of Magnetic Resonance* **158** (2002) 164–168.
This quantity can be optimised with respect to zero- and first-order
phase correction for automatic (unsupervised) phase correction.
"""
function entropy(s::SpectData{T,1}) where {T<:Number}
    h= s.dat .|> real  .|> abs
    h/=sum(h)
    return sum(ent.(h))
end



struct AutoPhaseCorrectChen <: NMRProcessor1D
    dim::Int64
    verbose::Bool
    γ::Float64
end

@doc raw"""
    function AutoPhaseCorrectChen(; dim::Integer=1, verbose::Bool=false, γ::Real=1.0e-5)

returns a processor that performs automatic phase correction of a spectrum along
the dimension `dim` (default `1`) using the minimum entropy algorithm by Chen et al. in
*Journal of Magnetic Resonance* **158** (2002) 164–168. The parameter `γ` can be
used to add a penalty term to the optimisation target, which penalises negative
peaks in the spectrum. This can be useful to avoid overcorrection in noisy
spectra.
"""
function AutoPhaseCorrectChen(; dim::Integer=1, verbose::Bool=false, γ::Real=1.0e-5)
    return AutoPhaseCorrectChen(dim,verbose,Float64(γ))
end

# penalty(x) computes the sum of squares of all negative points in x
function penalty(x)
    # x /= sum(abs.(x))
    return sum(k<0.0 ? k*k : 0.0 for k in x)
end

# this is the minimisation target for automatic phase correction
function goalfun(x,spect,γ)
    pc = PhaseCorrect(ph0=x[1], ph1=x[2]/1000)
    c = pc(spect) 
    return entropy(c)+γ*penalty(real.(c.dat));
end

function (apc::AutoPhaseCorrectChen)(spect::SpectData{T,1}) where {T<:Number}
    dspect = Derivative()(spect)
    # do a 1D optimisation of the zero-order pc first 
    res0  = Optim.optimize(x->goalfun([x[1],0.0],dspect,apc.γ),-pi,pi,Optim.Brent());
    
    if apc.verbose print(res0) end;        
    p0 = Optim.minimizer(res0)[1]

    # then a 2D optimisation of zero- and first-order pc
    result=Optim.optimize(x->goalfun(x,dspect,apc.γ),[p0,0.0],
            Optim.BFGS(),
            Optim.Options(show_trace=false,
                          f_calls_limit=500,
                          time_limit=2.0,
                          g_tol=1.0e-8)
        );
    if apc.verbose print(result) end;
    pc=PhaseCorrect(ph0=Optim.minimizer(result)[1], ph1=Optim.minimizer(result)[2]/1000);
    scorr = pc(spect);
  return scorr ;
end


@doc raw"""
    function PeakAlign(; dim::Integer=1, readpos::Real, wdw::Integer)

returns a processor that aligns a spectrum along the dimension `dim`
(default `1`) to a specific position `readpos`.  It works by finding a
maximum in the spectrum within a window `wdw` that
is closest to `readpos`, and then shifting the spectrum such that this maximum
is exactly at `readpos`. This can be useful to align spectra to a reference
peak, e.g., TMS.
"""
struct PeakAlign <: NMRProcessor1D
    dim::Int64
    readpos::Float64
    wdw::Int64
    PeakAlign(; dim::Integer=1, readpos::Real, wdw::Integer) = new(dim, Float64(readpos), wdw)
end

function (pa::PeakAlign)(spect::SpectData{T,1}) where {T<:Number}
    # find the index of the point closest to readpos
    idx = findmin(abs.(spect.coord[1] .- pa.readpos))[2]
    # find a maximum in the spectrum that is closest to readpos
    maxidx = findmax(abs.(spect.dat[idx-pa.wdw:idx+pa.wdw]))[2] + idx - pa.wdw
    newdat = circshift(spect.dat, idx-maxidx)
    return SpectData(newdat, (spect.coord[1],))
end

@doc raw"""
    function DigitalFilter(b::Vector{ComplexF64}; dim::Integer=1)

returns a processor that applies a digital FIR filter with coefficients `b`
along the dimension `dim` (default `1`), via `DSP.filt`. The coefficients
`b` can be produced by `BandReject` or `BandPass`.
"""
struct DigitalFilter <: NMRProcessor1D
    b::Vector{ComplexF64}
    # a::Vector{Float64}
    dim::Int64
    DigitalFilter(b::Vector{ComplexF64}; dim::Integer=1) = new(b, dim)
end

import DSP

function (df::DigitalFilter)(spect::SpectData{T,1}) where {T<:Number}
    newdat = DSP.filt(df.b, spect.dat)
    return SpectData(newdat, (spect.coord[1],))
end

@doc raw"""
    function BandReject(lf::Float64, hf::Float64, n::Integer)

returns the coeffiecients of a  digital band-rejection filter with lower and upper cutoff frequencies `lf`
and `hf`, respectively, and filter order `n`. 
The frequencies are given as a fraction of the spectral width. 
The filter is designed using the
window method, with a Blackman window. The returned filter coefficients can be
used to create a `DigitalFilter` processor. 
"""
function BandReject(lf,hf,n)
   b = [ t == 0 ? -ComplexF64(hf-lf,0.0) : 1.0/(2pi*im*t)*(exp(2pi*im*lf*t)-exp(2pi*im*hf*t)) for t=-n:n ]
   b .*= 0.42 .- 0.5*cos.(pi/n*(0:2n)) .+ 0.08*cos.(2pi/n*(0:2n))
   b[n+1] += ComplexF64(1.0,0.0)
   return b
end

@doc raw"""
    function BandPass(lf::Float64, hf::Float64, n::Integer)

returns the coeffiecients of a  digital band-pass filter with lower and upper cutoff frequencies `lf`
and `hf`, respectively, and filter order `n`. 
The frequencies are given as a fraction of the spectral width. 
The filter is designed using the
window method, with a Blackman window. The returned filter coefficients can be
used to create a `DigitalFilter` processor. 
"""
function BandPass(lf,hf,n)
   b = -[ t == 0 ? -ComplexF64(hf-lf,0.0) : 1.0/(2pi*im*t)*(exp(2pi*im*lf*t)-exp(2pi*im*hf*t)) for t=-n:n ]
   b .*= 0.42 .- 0.5*cos.(pi/n*(0:2n)) .+ 0.08*cos.(2pi/n*(0:2n))
   return b
end
