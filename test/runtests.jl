using NMRflux
using Test
using LinearAlgebra
using Statistics

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
bruker_file() = NMRflux.Examples.Data["HCC cell culture media spectra"]["files"][1]
jeol_file()   = NMRflux.Examples.Data["DMEM culture medium"]["files"][1]

# Full processing pipeline for Bruker FID -> spectrum
function bruker_spectrum(n=2^16)
    params, d = NMRflux.load(bruker_file(), :Bruker)
    zf  = ZeroFill([n])
    ft  = FourierTransform([n], [1])
    d |> zf |> ft, params, d
end

# ---------------------------------------------------------------------------
# 1. SpectData construction & accessors
# ---------------------------------------------------------------------------
@testset "SpectData construction and accessors" begin
    A = ones(40, 50)
    c = (range(0.0, 1.0, 40), range(1.0, 2.0, 50))
    S = NMRflux.SpectData(A, c)

    @test size(S) == (40, 50)
    @test coords(S) === S.coord
    @test coords(S, 1) === c[1]
    @test coords(S, 2) === c[2]

    # inferred coords
    S2 = NMRflux.SpectData(rand(10, 8))
    @test size(S2) == (10, 8)
    @test length(coords(S2, 1)) == 10
    @test length(coords(S2, 2)) == 8

    # scalar element access
    @test S[1, 1] == 1.0
    @test S[5, 7] == 1.0
end

# ---------------------------------------------------------------------------
# 2. SpectData slicing & broadcasting
# ---------------------------------------------------------------------------
@testset "SpectData slicing and broadcasting" begin
    A = reshape(1.0:60.0, 3, 4, 5)
    c = (range(0.0,1.0,3), range(0.0,1.0,4), range(0.0,1.0,5))
    S = NMRflux.SpectData(collect(A), c)

    # slice out one index in first dim -> 2D SpectData
    sl = S[1, :, :]
    @test isa(sl, NMRflux.SpectData{Float64,2})
    @test size(sl) == (4, 5)
    @test coords(sl, 1) === c[2]
    @test coords(sl, 2) === c[3]

    # slice out one index in second dim -> 2D SpectData
    sl2 = S[:, 1, :]
    @test isa(sl2, NMRflux.SpectData{Float64,2})
    @test size(sl2) == (3, 5)

    # broadcasting: .* scalar preserves SpectData type and coords
    S1d = NMRflux.SpectData(collect(1.0:10.0), (range(0.0,1.0,10),))
    R = S1d .* 2.0
    @test isa(R, NMRflux.SpectData)
    @test R.coord[1] === S1d.coord[1]
    @test all(R.dat .≈ S1d.dat .* 2.0)

    # broadcasting: S .+ S
    R2 = S1d .+ S1d
    @test isa(R2, NMRflux.SpectData)
    @test all(R2.dat .≈ S1d.dat .* 2.0)
end

@testset "Craft.analyze" begin
    ex = NMRflux.Examples.Data["HCC cell culture media spectra"]["files"][1]
    params, fid_full = NMRflux.load(ex, :Bruker)
    fid = fid_full[1:4096]   # short section, for a fast test

    K = 20
    result = NMRflux.Craft.analyze(fid, K)
    tbl = result["resonance_table"]

    SW = 1.0 / step(NMRflux.coords(fid, 1))

    @test size(tbl, 1) == K
    @test all(-SW/2 .<= tbl.frequency_hz .<= SW/2)
    @test all(isfinite, tbl.intensity)
    @test all(isfinite, tbl.phase_rad)
    @test all(isfinite, tbl.linewidth_hz)
end

# ---------------------------------------------------------------------------
# 3. FileIO — Bruker
# ---------------------------------------------------------------------------
@testset "FileIO Bruker" begin
    params, d = NMRflux.load(bruker_file(), :Bruker)
    @test isa(params, Dict)
    @test haskey(params, "SW_h")
    @test haskey(params, "GRPDLY")

    @test isa(d, NMRflux.SpectData{ComplexF64,1})
    @test size(d) == (length(coords(d,1)),)
    @test length(d) > 0

    # regression checksum (computed from first run, kept as guard)
    @test isapprox(sum(d), 2.616357072520625e11 + 1.335611666755e11im; rtol=1e-6)

    # unsupported vendor must throw
    @test_throws Exception NMRflux.load(bruker_file(), :UnknownVendor)
end

# ---------------------------------------------------------------------------
# 4. FileIO — JEOL 1D
# ---------------------------------------------------------------------------
@testset "FileIO JEOL 1D" begin
    params, d = NMRflux.load(jeol_file(), :JEOL)

    @test isa(params, Dict)
    @test haskey(params, "X_SWEEP")

    @test isa(d, NMRflux.SpectData{ComplexF64,1})
    @test size(d) == (length(coords(d,1)),)
    @test length(d) > 0

    # data should be genuinely complex
    @test norm(imag.(d.dat)) > 0.0
end

# ---------------------------------------------------------------------------
# 4b. FileIO — JEOL 2D (SH3 HSQC)
# ---------------------------------------------------------------------------
@testset "FileIO JEOL 2D" begin
    # The HSQC is the first file in the SH3 example set
    hsqc_file = NMRflux.Examples.Data["SH3 domain HSQC"]["files"][1]
    params, d = NMRflux.load(hsqc_file, :JEOL)

    # Should be a 2D SpectData
    @test isa(d, NMRflux.SpectData{ComplexF64,2})

    # Shape: 1280 direct (¹H) × 256 indirect (¹⁵N)
    @test size(d) == (1280, 256)
    @test size(d) == (length(coords(d,1)), length(coords(d,2)))

    # Both axes are time-domain: start at 0, end > 0
    @test first(coords(d,1)) == 0.0
    @test last(coords(d,1))  >  0.0
    @test first(coords(d,2)) == 0.0
    @test last(coords(d,2))  >  0.0

    # Direct axis (¹H at 600 MHz) acquisition time ≈ 0.170 s
    @test isapprox(last(coords(d,1)), 0.17046512; rtol=1e-4)

    # Data is genuinely complex in both real and imaginary parts
    @test norm(real.(d.dat)) > 0.0
    @test norm(imag.(d.dat)) > 0.0
end

# ---------------------------------------------------------------------------
# 5. ZeroFill
# ---------------------------------------------------------------------------
@testset "ZeroFill" begin
    params, d = NMRflux.load(bruker_file(), :Bruker)
    n_orig = length(d)
    n_new  = 2^16
    zf = ZeroFill([n_new])
    zd = zf(d)

    @test length(zd) == n_new
    @test length(coords(zd,1)) == n_new

    # first n_orig points are preserved
    @test all(zd.dat[1:n_orig] .== d.dat[1:n_orig])

    # remainder is zero
    @test all(zd.dat[n_orig+1:end] .== 0)
end

# ---------------------------------------------------------------------------
# 6. FourierTransform
# ---------------------------------------------------------------------------
@testset "FourierTransform" begin
    params, d = NMRflux.load(bruker_file(), :Bruker)
    n = 2^16
    zf = ZeroFill([n])
    ft = FourierTransform([n], [1])
    sp = ft(zf(d))

    @test length(sp) == n
    @test length(coords(sp,1)) == n

    # freq axis is symmetric around 0
    @test first(coords(sp,1)) < 0.0
    @test last(coords(sp,1))  > 0.0

    @test eltype(sp.dat) == ComplexF64

    # fftshift=false: coord is still the symmetric range (same as shifted),
    # but the data array itself differs because no circular shift was applied
    ft_no_shift = FourierTransform([n], [1]; fftshift=false)
    sp2 = ft_no_shift(zf(d))
    @test length(sp2) == n
    @test coords(sp2,1) == coords(sp,1)          # coord axis unchanged
    @test !all(sp2.dat .≈ sp.dat)                # data differs (no fftshift applied)
end

# ---------------------------------------------------------------------------
# 7. Apodize
# ---------------------------------------------------------------------------
@testset "Apodize" begin
    params, d = NMRflux.load(bruker_file(), :Bruker)

    # R=0 is identity
    ap0 = Apodize([0.0])
    d0  = ap0(d)
    @test all(d0.dat .≈ d.dat)
    @test coords(d0,1) === coords(d,1)

    # large R strongly attenuates tail relative to head
    ap_big = Apodize([50.0])
    d_big  = ap_big(d)
    ratio_first = abs(d_big.dat[1]) / max(abs(d.dat[1]), 1e-30)
    ratio_last  = abs(d_big.dat[end]) / max(abs(d.dat[end]), 1e-30)
    @test ratio_last < ratio_first   # tail suppressed more than head
end

# ---------------------------------------------------------------------------
# 8. PhaseCorrect
# ---------------------------------------------------------------------------
@testset "PhaseCorrect" begin
    sp, params, _ = bruker_spectrum()

    # identity phase (ph0=0, ph1=0)
    pc0 = PhaseCorrect(0.0, 0.0, 1)
    sp0 = pc0(sp)
    @test all(sp0.dat .≈ sp.dat)
    @test coords(sp0,1) === coords(sp,1)

    # ph0=pi negates the real part of a purely real synthetic signal
    v = NMRflux.SpectData(ones(ComplexF64, 100), (range(0.0,1.0,100),))
    pc_pi = PhaseCorrect(pi, 0.0, 1)
    vr = pc_pi(v)
    @test all(real.(vr.dat) .≈ -1.0)
end

# ---------------------------------------------------------------------------
# 9. Derivative
# ---------------------------------------------------------------------------
@testset "Derivative" begin
    # derivative of a constant synthetic signal is (approximately) zero
    v = NMRflux.SpectData(ones(ComplexF64, 200), (range(0.0, 1.0, 200),))
    der = Derivative(1)
    dv  = der(v)
    @test length(dv) == length(v)
    @test coords(dv,1) === coords(v,1)
    # interior points should be essentially zero
    @test maximum(abs.(dv.dat[5:end-4])) < 1e-10

    # derivative of a unit ramp is approximately 1.0 everywhere (interior)
    ramp = NMRflux.SpectData(collect(ComplexF64, range(0.0, 1.0, 200)),
                              (range(0.0, 1.0, 200),))
    dr = der(ramp)
    interior = real.(dr.dat[5:end-4])
    @test all(abs.(interior .- 1.0) .< 0.01)
end

# ---------------------------------------------------------------------------
# 10. Integral
# ---------------------------------------------------------------------------
@testset "Integral" begin
    # integral of a constant = linearly ramped cumsum
    n   = 100
    inc = 1.0/(n-1)
    v   = NMRflux.SpectData(ones(ComplexF64, n), (range(0.0, 1.0, n),))
    int = Integral(1)
    iv  = int(v)

    @test length(iv) == n
    @test coords(iv,1) === coords(v,1)

    expected = cumsum(ones(n)) .* inc
    @test all(real.(iv.dat) .≈ expected)

    # Integral of a sinusoid over one full period is approximately zero
    sine = NMRflux.SpectData(sin.(range(0.0, 2pi, 200)) .+ 0im,
                              (range(0.0, 2pi, 200),))
    iv2  = int(sine)
    @test abs(real(iv2.dat[end])) < 0.1   # integral of sin over [0,2π] ≈ 0
end

# ---------------------------------------------------------------------------
# 11. MedianBaselineCorrect
# ---------------------------------------------------------------------------
@testset "MedianBaselineCorrect" begin
    params, d = NMRflux.load(bruker_file(), :Bruker)
    n = 2^16
    zf  = ZeroFill([n])
    ft  = FourierTransform([n],[1])
    pc  = PhaseCorrect(0.80pi, 2pi*0.00172, 1)
    bc  = MedianBaselineCorrect(1, wdw=1024)

    sp  = d |> zf |> ft |> pc |> bc

    @test isa(sp, NMRflux.SpectData)
    @test length(sp) == n
    @test eltype(sp.dat) <: Real   # MedianBaselineCorrect returns real
end

# ---------------------------------------------------------------------------
# 12. AutoPhaseCorrectChen
# ---------------------------------------------------------------------------
@testset "AutoPhaseCorrectChen" begin
    params, d = NMRflux.load(bruker_file(), :Bruker)
    n   = 2^15   # must be >= raw FID length (~32693)
    zf  = ZeroFill([n])
    ft  = FourierTransform([n],[1])
    apc = AutoPhaseCorrectChen(1; verbose=false)

    sp  = d |> zf |> ft |> apc

    @test isa(sp, NMRflux.SpectData)
    @test length(sp) == n
    @test norm(real.(sp.dat)) > 0.0
end

# ---------------------------------------------------------------------------
# 13. PeakAlign
# ---------------------------------------------------------------------------
@testset "PeakAlign" begin
    # PeakAlign searches within `wdw` index steps of readpos for the maximum,
    # then circularly shifts the data so that maximum lands on readpos.
    # Peak must therefore start within wdw steps of readpos.
    n        = 1000
    crd      = range(0.0, 10.0, n)
    step_val = step(crd)
    peak_pos = 5.1     # slightly offset from readpos
    readpos  = 5.0
    wdw      = 50      # 50 * step ≈ 0.5 units — comfortably covers the 0.1 offset
    vals = ComplexF64[ 1.0 / (1.0 + ((x - peak_pos)/0.05)^2) for x in crd ]
    sp   = NMRflux.SpectData(vals, (crd,))

    pa  = PeakAlign(1, readpos, wdw)
    spa = pa(sp)

    aligned_max_pos = coords(spa,1)[argmax(abs.(spa.dat))]
    @test abs(aligned_max_pos - readpos) <= 2*step_val
end

# ---------------------------------------------------------------------------
# 14. Chain / pipeline composition
# ---------------------------------------------------------------------------
@testset "Chain and pipeline" begin
    params, d = NMRflux.load(bruker_file(), :Bruker)
    n = 2^15

    zf = ZeroFill([n])
    ft = FourierTransform([n],[1])

    # |> piping
    sp_pipe  = d |> zf |> ft

    # Chain (applies left-to-right because Chain reverses ∘)
    sp_chain = Chain(zf, ft)(d)

    @test length(sp_pipe)  == n
    @test length(sp_chain) == n
    @test all(sp_pipe.dat .≈ sp_chain.dat)

    # full end-to-end pipeline on Bruker data
    pc = PhaseCorrect(0.80pi, 2pi*0.00172, 1)
    bc = MedianBaselineCorrect(1, wdw=512)
    sp_full = d |> zf |> ft |> pc |> bc

    @test isa(sp_full, NMRflux.SpectData)
    @test length(sp_full) == n
end
