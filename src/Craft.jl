module Craft

using NMRflux
using LinearAlgebra
using LsqFit
using DataFrames

export analyze


function esprit(x, K)
    N = length(x)
    L = N ÷ 2
    M = N - L + 1

    X    = [x[i+j-1] for i in 1:M, j in 1:L]
    svdX = svd(X)
    Us   = svdX.U[:, 1:K]

    Phi = Us[1:M-1, :] \ Us[2:M, :]
    λ   = eigvals(Phi)

    ωs     = angle.(λ)       # rad/sample
    decays = -log.(abs.(λ))  # 1/sample, positive for damped sinusoids

    idx = sortperm(ωs)
    return ωs[idx], decays[idx], svdX.S
end

@doc raw"""
    esprit(data::SpectData{T,1}, K) where T

Apply the ESPRIT algorithm to 1D spectral data to estimate the frequencies and decay rates
of the underlying sinusoids.

The function returns the frequencies in Hz, decay rates in Hz, and the singular values of the Hankel matrix.

**Example:**
```julia
using NMRflux
fidData = randn(1024) .+ im*randn(1024)
data = SpectData(fidData, (0:1023,))
freqs, decays, S = esprit(data, 5)
```

ESPRIT (Estimation of Signal Parameters via Rotational Invariance Techniques) is a
subspace method for estimating the frequencies of damped sinusoids in noise. It exploits
the algebraic structure of the signal subspace rather than searching over a frequency grid,
giving superresolution estimates that are not limited by the FFT bin spacing.

**Algorithm outline:**

1. **Hankel matrix.** From the FID $x[1],\ldots,x[N]$, form the $M \times L$ Hankel matrix
   ($L = \lfloor N/2 \rfloor$, $M = N - L + 1$):
   $$X_{ij} = x[i+j-1]$$
   For $K$ ideal damped sinusoids, $X$ has rank $K$ in the absence of noise.

2. **Signal subspace.** Compute the SVD of $X$. The $K$ dominant left singular vectors
   span the signal subspace $U_s$ ($M \times K$).

3. **Rotational invariance.** Split $U_s$ into two overlapping sub-matrices differing by
   one row:
   $$U_s^\uparrow = U_s[1{:}M-1,\,:\,], \qquad U_s^\downarrow = U_s[2{:}M,\,:\,]$$
   Because consecutive rows of $X$ are related by a time shift, these satisfy
   $U_s^\downarrow \approx U_s^\uparrow \Phi$, where $\Phi$ is a $K\times K$ matrix
   whose eigenvalues are $e^{(-\delta_k + i\omega_k)\Delta t}$.

4. **Frequency extraction.** Solve for $\Phi$ by least squares
   ($\Phi = (U_s^\uparrow)^+ U_s^\downarrow$) and take the angles of its eigenvalues
   to recover the frequencies $\omega_k$ in rad/sample.

**Reference:**
R. Roy and T. Kailath, "ESPRIT — Estimation of Signal Parameters via Rotational Invariance
Techniques," *IEEE Trans. Acoust., Speech, Signal Process.*, vol. 37, no. 7,
pp. 984–995, July 1989.
"""
function esprit(data::SpectData{T,1}, K) where T
    freqs, decays, S = esprit(data.dat, K)
    dt       = step(NMRflux.coords(data, 1))
    freq_hz  = freqs ./ dt ./ (2π)
    decay_hz = decays ./ dt
    return freq_hz, decay_hz, S
end


@doc raw"""
      laplace_log_evidence(fit::LsqFit.LsqFitResult)

Compute the Laplace approximation to the log-evidence of a nonlinear
least-squares fit.

To compare models with different numbers of oscillators $K$, we need the
**model evidence** $p(\mathbf{x} \mid K)$ — the probability of the data after
marginalising over all parameter values:
$$p(\mathbf{x} \mid K) = \int p(\mathbf{x} \mid \mathbf{B}, K)\, p(\mathbf{B} \mid K)\, d\mathbf{B}$$
A model with more parameters fits better, but spreads its prior over a larger space,
incurring an automatic Occam penalty. The evidence balances fit quality against
model complexity without manual regularisation.

**Laplace approximation.** For a Gaussian likelihood with noise variance $\sigma^2$
and a broad, flat prior, expanding the log-posterior to second order around the
maximum-likelihood estimate $\hat{\mathbf{B}}$ gives:
$$\ln p(\mathbf{x} \mid K) \approx
  -\frac{n}{2}\ln\sigma^2
  +\frac{1}{2}\ln\det(J^\top J)
  -\frac{K}{2}\ln n$$
where $J$ is the Jacobian of the residuals evaluated at $\hat{\mathbf{B}}$,
$n$ is the number of data points, and $K$ is the number of fitted parameters.
The three terms are respectively the log-likelihood at the optimum, the log-volume
of the posterior peak (how well the data constrain the parameters), and the BIC
penalty for model complexity.

To select the number of oscillators, fit for $K = 1, 2, \ldots$ and choose the
$K$ with the highest log evidence.

**Numerical note.** When some parameters are poorly constrained (e.g. duplicate
frequencies), $J^\top J$ becomes singular. A small ridge $\varepsilon I$ is added
before taking the log-determinant to regularise this case.
"""
function laplace_log_evidence(fit)
    K   = length(fit.param)
    n   = length(fit.resid)
    σ²  = sum(fit.resid.^2) / n
    JtJ = fit.jacobian' * fit.jacobian
    ε   = 1e-10 * maximum(diag(JtJ))
    return -n/2 * log(σ²) + 1/2 * logabsdet(JtJ + ε*I)[1] - K/2 * log(n)
end

@doc raw"""
      oscillators_fixed_freq(r, B, freq_hz)
Compute the sum of $K$ damped complex exponentials with fixed frequencies.

The FID is modelled as a superposition of $K$ damped complex sinusoids:
$$x(t) = \sum_{k=1}^{K} c_k \, e^{(-\delta_k + i\omega_k)t}$$
where $c_k = a_k e^{i\phi_k}$ is a complex amplitude encoding both the peak
amplitude $a_k$ and phase $\phi_k$, $\omega_k = 2\pi f_k$ is the angular
frequency, and $\delta_k > 0$ is the decay rate (inversely proportional to the
transverse relaxation time $T_2$).

**Why fix the frequencies?**
In principle one could fit all $4K$ parameters $(\text{Re}\,c_k,\,\text{Im}\,c_k,\,\omega_k,\,\delta_k)$
simultaneously. In practice this leads to a severely ill-conditioned optimisation
for two reasons.
First, the cost function has many nearly equivalent local minima: two oscillators
at slightly different frequencies can produce very similar residuals, so the
optimiser can swap, merge, or drift them without penalty.
Second, the Jacobian columns for $\omega_k$ and $\delta_k$ are proportional to
$c_k \cdot t \cdot e^{(-\delta_k+i\omega_k)t}$; if two oscillators converge to
the same frequency their columns become linearly dependent and the Jacobian
becomes singular.

ESPRIT provides accurate, superresolution frequency estimates directly from the
data, so there is no reason to re-estimate them by nonlinear optimisation.
Fixing $\omega_k$ at the ESPRIT values reduces the problem to $3K$ parameters
$( \text{Re}\,c_k,\,\text{Im}\,c_k,\,\delta_k )$.
For fixed frequencies the model is **linear** in the complex amplitudes $c_k$,
and only mildly nonlinear through the decay rates $\delta_k$, giving a
much better conditioned fit.

The decay rate is reparametrised as $\delta_k = \delta_{\max}\,\sigma(\gamma_k)$
where $\sigma$ is the sigmoid function, bounding $\delta_k$ strictly within
$(0, \delta_{\max})$ and preventing runaway decay during optimisation.
"""
function oscillators_fixed_freq(r, B, freq_hz)
    K        = length(freq_hz)
    freq_rads = freq_hz .* 2π
    return sum(1:K) do k
        c     = B[3k-2] + im*B[3k-1]
        decay = decay_rate_from_param(B[3k])
        c .* exp.((-decay + im*freq_rads[k]) .* r)
    end
end

const DECAY_MAX = 40π   # ~20 Hz maximum linewidth
function decay_rate_from_param(γ)
    return DECAY_MAX / (1 + exp(-γ))
end
function param_from_decay_rate(δ)
    return log(δ / (DECAY_MAX - δ))
end


@doc raw"""
      oscillators_fixed_freq_jacobian(r, B, freq_hz)

Compute the Jacobian of the sum of $K$ damped complex exponentials with fixed
frequencies.
cf. `oscillators_fixed_freq`.
"""
function oscillators_fixed_freq_jacobian(r, B, freq_hz)
    K         = length(freq_hz)
    freq_rads = freq_hz .* 2π
    nr        = length(r)
    J         = zeros(2nr, 3K)
    for k in 1:K
        c     = B[3k-2] + im*B[3k-1]
        decay = decay_rate_from_param(B[3k])
        e     = exp.((-decay + im*freq_rads[k]) .* r)
        ce    = c .* e

        J[1:nr,     3k-2] =  real(e);  J[nr+1:end, 3k-2] = imag(e)
        J[1:nr,     3k-1] = -imag(e);  J[nr+1:end, 3k-1] = real(e)

        d_factor = decay * (1 - decay / DECAY_MAX)
        d = -r .* ce .* d_factor
        J[1:nr,     3k]   =  real(d);  J[nr+1:end, 3k]   = imag(d)
    end
    return J
end

@doc raw"""
      fit_fid_fixed_freq(fid, r, freq_hz, decay_rate)

Fit a sum of $K$ damped complex exponentials with fixed frequencies to an FID
`fid` sampled at times `r`. The frequencies are fixed at `freq_hz`, and the
decay rates are initialised at `decay_rate` (in Hz).

Returns a `Dict` with keys:
- `"n_oscillators"`   — $K$
- `"residual"`        — sum of squared residuals of the fit
- `"log_evidence"`    — Laplace approximation to the log-evidence (see `laplace_log_evidence`)
- `"fit"`             — the underlying `LsqFit.LsqFitResult`, giving access to
                         fit diagnostics such as confidence intervals
- `"resonance_table"` — a `DataFrame`, sorted by frequency, with columns
                         `frequency_hz`, `intensity`, `phase_rad` and `linewidth_hz`
- `"resonances"`      — the same resonances as a vector of `Dict`s, each with keys
                         `frequency_hz`, `amplitude`, `phase_rad`, `decay_rate_hz` and `T2_s`

cf. `oscillators_fixed_freq` and `oscillators_fixed_freq_jacobian`.
"""
function fit_fid_fixed_freq(fid::AbstractVector{<:Complex}, r::StepRangeLen,
                             freq_hz::Vector, decay_rate::Vector)
    K  = length(freq_hz)
    dt = step(r)

    model(r, B)    = let pred = oscillators_fixed_freq(r, B, freq_hz)
        [real(pred .- fid); imag(pred .- fid)]
    end
    jacobian(r, B) = oscillators_fixed_freq_jacobian(r, B, freq_hz)

    # Initialise γ from ESPRIT decay estimates via sigmoid inverse (logit)
    decay_clamped = clamp.(decay_rate, 1e-3, DECAY_MAX - 1e-3)
    γ0 = log.(decay_clamped ./ (DECAY_MAX .- decay_clamped))

    B0          = zeros(3K)
    B0[1:3:end] .= 1.0   # initial cr
    B0[3:3:end]  = γ0

    fit = curve_fit(model, jacobian, r, zeros(2length(r)), B0)

    B          = fit.param
    amplitude  = abs.(B[1:3:end] .+ im .* B[2:3:end])
    phase_rad  = angle.(B[1:3:end] .+ im .* B[2:3:end])
    decay_rate_fit = DECAY_MAX ./ (1 .+ exp.(-B[3:3:end]))

    idx = sortperm(freq_hz)
    result = Dict(
        "n_oscillators" => K,
        "residual"      => sum(fit.resid.^2),
        "log_evidence"  => laplace_log_evidence(fit),
        "fit"            => fit,
        "resonance_table" => DataFrame(
            frequency_hz = freq_hz[idx],
            intensity    = amplitude[idx],
            phase_rad    = phase_rad[idx],
            linewidth_hz = decay_rate_fit[idx] ./ (2π),
        ),
        "resonances"    => [
            Dict(
                "frequency_hz"  => freq_hz[i],
                "amplitude"     => amplitude[i],
                "phase_rad"     => phase_rad[i],
                "decay_rate_hz" => decay_rate_fit[i] / (2π),
                "T2_s"          => π / decay_rate_fit[i]
            )
            for i in idx
        ]
    )

    return result
end

@doc raw"""
      fit_fid_fixed_freq(fid::SpectData{T,1}, freq_hz::Vector, decay_rate::Vector) where T

Fit a sum of $K$ damped complex exponentials with fixed frequencies to an FID
`fid` sampled at times `r`. The frequencies are fixed at `freq_hz`, and the
decay rates are initialised at `decay_rate` (in Hz).

Returns the same `Dict` as `fit_fid_fixed_freq(fid::AbstractVector{<:Complex}, r, freq_hz, decay_rate)`
— see that method's docstring for the list of keys.

cf. `oscillators_fixed_freq` and `oscillators_fixed_freq_jacobian`.
"""
function fit_fid_fixed_freq(fid::SpectData{T,1}, freq_hz::Vector, decay_rate::Vector) where T
    result = fit_fid_fixed_freq(fid.dat, NMRflux.coords(fid, 1), freq_hz, decay_rate)
    return result
end

@doc raw"""
      analyze(fid::SpectData{T,1}, K::Integer=150) where T

Decompose an FID into a sum of damped sinusoids: estimate `K` resonance
frequencies and decay rates via `esprit`, and refine them with
`fit_fid_fixed_freq`.

Returns the `Dict` produced by `fit_fid_fixed_freq` — see its docstring for
the full list of keys. The resonance table and fit diagnostics of most
interest are found under the `"resonance_table"` and `"fit"` keys, respectively.
"""
function analyze(fid::SpectData{T,1}, K::Integer=150) where T
    freq_hz, decays, _ = esprit(fid, K)
    result = fit_fid_fixed_freq(fid, freq_hz, decays)
    return result 
end


end
