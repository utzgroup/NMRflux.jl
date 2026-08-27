# RINSE

RINSE is a machine learning restoration demonstrator built on top of `NMRflux.jl`. It is not part of the NMRflux package. RINSE is developed and distributed separately, with its own repository, training code, inference code, and data generation utility. The purpose of this page is to show that the public simulation, data representation, and processing interfaces carry a complete external learning workflow from end to end. Nothing described here is part of the NMRflux API. The method itself, its generative model, and its results are reported in the RINSE publication and its supporting information.

!!! note "Dataset availability"
    The synthetic datasets used to train and evaluate the RINSE demonstrator are publicly available on Zenodo:

    **Kaleel, M. and Utz, M. (2026).**  
    *RINSE simulated 1H NMR spectra dataset (v1.0.0).*  
    Zenodo. [https://doi.org/10.5281/zenodo.20701118](https://doi.org/10.5281/zenodo.20701118)

    Please cite the Zenodo record if you use these datasets.

## 1. Scope

RINSE targets systematic artefacts in one dimensional NMR spectra: phase distortion, baseline corruption, and residual solvent or carrier contamination. It learns a mapping from an artefact degraded spectrum to the corresponding reference spectrum. The restored spectra are intended to improve interpretability and preprocessing quality. The RINSE demonstrator shows that an NMR specific machine learning workflow can be built on the NMRflux interfaces without leaving Julia and without an intermediate data format.

## 2. What RINSE takes from NMRflux

`NMRflux.SpinSim` provides the spin physics behind the synthetic training signals. `SpectData` is the container that moves data between simulation, processing, training, and inference, so no stage needs a conversion step or a bespoke array layout. The classical processing chain turns simulated free induction decays into the spectra the network sees, using the same objects documented elsewhere in this manual:

```julia
process = Chain(
    ZeroFill([2^16]),
    Apodize([0.5pi]),
    FourierTransform([2^16], [1]; fftshift=true)
)

spectrum = fid |> process
```

The first time domain point is halved before transformation, following the usual half point correction for discrete Fourier processing of FIDs. Every spectrum entering the learning pipeline therefore shares a length, a frequency grid, and a processing history. Everything else belongs to RINSE. The paired synthetic signals are produced by an application specific generator that ships with the RINSE repository, calls `SpinSim` for the underlying spin dynamics, and is **not part of the NMRflux core API**. Its generative model, its artefact models, and the composition of the released datasets are described in the supporting information of the RINSE publication, which is the authoritative source for them.

## 3. Learning formulation

Restoration is posed as supervised regression on pairs of spectra derived from the same underlying simulated signal. Because both members of a pair come from identical spin system parameters and pass through the same processing chain, the difference between them is the applied perturbation and nothing else. NMR spectra are complex valued, but the learning problem is expressed with real valued tensors so that standard layers apply. The network reads two channels, the real and the imaginary part of the degraded spectrum, and writes a single channel, the absorptive real part of the reference spectrum. In tensor form, with `W` spectral points and `N` examples in a batch,

```math
\mathbf{x} \in \mathbb{R}^{W \times 2 \times N}
\;\longmapsto\;
\hat{\mathbf{y}} \in \mathbb{R}^{W \times 1 \times N}.
```

Keeping the imaginary part on the input side preserves the phase information carried by the signal without requiring complex valued layers. Restricting the output to the real channel reflects where one dimensional NMR analysis is actually done. Every loss term and every reported metric is computed on that real valued output.

### 3.1 Normalisation

Spectral amplitudes vary by orders of magnitude across concentration, receiver gain, and scaling convention. Passing that range straight to a network lets the largest signals dominate the gradient, so every spectrum is normalised before it reaches the model. A single scalar is computed from the degraded complex spectrum,

```math
\mathcal{N}(\hat{\mathbf{S}}) = \max_k \left| \hat{S}_k \right|,
```

taken over the magnitudes of the complex points before any split into channels. That same scalar normalises both members of the pair,

```math
\tilde{\mathbf{S}}_{\mathrm{degraded}}
=
\frac{\hat{\mathbf{S}}_{\mathrm{degraded}}}{\mathcal{N}(\hat{\mathbf{S}}_{\mathrm{degraded}})},
\qquad
\tilde{\mathbf{S}}_{\mathrm{reference}}
=
\frac{\hat{\mathbf{S}}_{\mathrm{reference}}}{\mathcal{N}(\hat{\mathbf{S}}_{\mathrm{degraded}})}.
```

Both members share one factor so that the amplitude relationship between them survives. The degraded input acquires unit maximum complex magnitude, and the reference keeps its amplitude relative to that input instead of being rescaled independently. Prediction then happens in the normalised domain and the result is returned to the original scale afterwards,

```math
\tilde{\mathbf{y}} = \mathcal{M}(\tilde{\mathbf{S}}_{\mathrm{degraded}}),
\qquad
\hat{\mathbf{y}} = \mathcal{N}(\hat{\mathbf{S}}_{\mathrm{degraded}})\,\tilde{\mathbf{y}},
```

where the calligraphic `M` denotes the trained network. Alternatives to the maximum magnitude scheme exist. An l2 norm is less brittle against an isolated spike but depends on the spectral grid, so zero filling or a change of length alters it, and a single reference bin is simple but fails when that bin sits near zero.

## 4. Predictor architecture

RINSE uses a one dimensional U-Net designed for paired spectral restoration. The encoder has three resolution levels followed by a bottleneck, and the decoder mirrors it. At each encoder level the input passes through a convolution of kernel width 3, then batch normalisation and a rectified linear unit, and then a residual refinement block of two further convolutional layers with batch normalisation and a skip connection. Writing the block input as `z`, its output is

```math
\mathbf{z}_{\mathrm{out}} = \mathbf{z} + \mathcal{F}(\mathbf{z}),
```

with `F` the learned residual mapping. Local spectral corrections accumulate through those blocks without the network having to relearn the identity at every level, which keeps optimisation stable. Downsampling between stages is mean pooling with factor 2, and the channel width grows as the spectral width shrinks: 32 channels at the first level, then 64, then 128, and 256 in the bottleneck. The bottleneck itself is a 3 point convolution followed by two residual blocks, and it holds the widest receptive field in the network. Each decoder stage upsamples by nearest neighbour interpolation with factor 2, applies a 3 point convolution with batch normalisation and a rectified linear unit, concatenates the corresponding encoder feature map through a skip connection, and passes the result through a block of the same form as the encoder. Those skip connections do the work that matters for NMR. They carry fine local detail past the pooling stages, so narrow resonances and sharp peak boundaries survive while the deeper layers still see broad spectral context. A final 1 by 1 convolution maps the last decoder representation to the single output channel. Since the three downsampling stages require a length divisible by 8, the spectral width is padded when necessary and cropped back after the output convolution, which leaves the architecture applicable to any one dimensional length while preserving exact output size. The network holds approximately 1.44 million trainable parameters.

## 5. Composite training objective

A plain pointwise regression loss concentrates optimisation pressure on dominant peaks and broad high amplitude regions, which leaves it insensitive to weak but chemically meaningful features. The objective used for the reported model therefore combines a reconstruction term with three NMR informed constraints,

```math
\mathcal{L}_{\mathrm{total}}
=
\mathcal{L}_{\mathrm{MSE}}
+
\lambda_{\mathrm{peak}}\,\mathcal{L}_{\mathrm{peak}}
+
\lambda_{\mathrm{L1}}\,\mathcal{L}_{\mathrm{L1}}
+
\lambda_{\mathrm{integ}}\,\mathcal{L}_{\mathrm{integ}},
```

where the coefficients set the contribution of the three auxiliary terms. The mean squared error anchors the optimisation and holds overall numerical agreement between prediction and reference,

```math
\mathcal{L}_{\mathrm{MSE}}
=
\frac{1}{M}\sum_{i=1}^{M}(\hat{y}_i-y_i)^2,
```

with `M` the number of scalar elements in the batch. A curvature weighted mean squared error then biases the optimisation toward narrow peak like structure,

```math
\mathcal{L}_{\mathrm{peak}}
=
\frac{1}{M}\sum_{i=1}^{M} w_i\,(\hat{y}_i-y_i)^2.
```

The weights come from the local negative curvature of the reference spectrum, estimated by the discrete second derivative

```math
d^{(2)}_i = y_{i-1} - 2y_i + y_{i+1},
```

so that

```math
w_i = 1 + \alpha\,c_i,
```

where `c` is the positive part of the negated second derivative, normalised per spectrum to lie between 0 and 1, and `alpha` controls the strength of the emphasis. Points at the top of a resonance carry strong negative curvature and are weighted up. Flat regions are not. An L1 residual term penalises deviations linearly instead of quadratically, which reduces the influence of isolated large residuals on a sparse spectrum,

```math
\mathcal{L}_{\mathrm{L1}}
=
\frac{1}{M}\sum_{i=1}^{M} \left| \hat{y}_i-y_i \right|.
```

Downstream NMR analysis depends on integrated signal, not only on pointwise agreement, so the spectrum is finally partitioned into `K` contiguous bins and the predicted and reference integrals are compared within each,

```math
\mathcal{L}_{\mathrm{integ}}
=
\frac{1}{K}\sum_{k=1}^{K}
\left(
\frac{\hat{I}_k-I_k}{S}
\right)^2,
\qquad
\hat{I}_k = \sum_{i \in \Omega_k} \hat{y}_i,
\qquad
I_k = \sum_{i \in \Omega_k} y_i,
```

where the interval of bin `k` is written as a subscripted omega, and `S` is the mean absolute amplitude of the reference spectrum multiplied by its length. That quantity is close to the total absolute integral, so the ratio is dimensionless and the term stays comparable across spectra of different amplitude. The coefficients in section 6 set the balance between the terms. `lambda_peak` is the largest at 5.0, which puts most of the auxiliary weight on resonance structure, and `lambda_l1` is the smallest at 0.05.

## 6. Training configuration

The coefficients above and the optimiser settings are read at run time from a TOML file, so the reported run reproduces without editing source code.

```toml
[Trainer]
    batchSize      = 64
    epochs         = 500
    timestamps     = 10
    eta            = 1.0e-4
    etaDecay       = 0.0
    beta1          = 0.9
    beta2          = 0.999
    weightDecay    = 1.0e-5
    fname          = "unet_nmr"

[Loss]
    name         = "nmr_peak_loss"
    alpha        = 10.0
    lambda_peak  = 5.0
    lambda_l1    = 0.05
    lambda_integ = 0.5
    K            = 32
    eps          = 1.0e-6
```

`[Trainer]` sets the AdamW hyperparameters and the checkpoint interval. `breakCriterion` is carried in the configuration file but the current training loop does not read it, so a run continues to the epoch budget unless it is stopped by hand. `[Loss]` sets the coefficients of the composite objective: `alpha` scales the curvature weighting, the three `lambda` entries weight the auxiliary terms, `K` is the number of contiguous bins used by the integral term, and `eps` is the floor added to the scale factor in the integral term. The learning rate follows an inverse time schedule,

```math
\eta_t = \frac{\eta}{1 + t\,\lambda_{\mathrm{decay}}}
```

at epoch `t`, where the decay constant is the `etaDecay` entry above. With `etaDecay` at zero the rate stays constant for the reported run, which was stopped by hand at 400 epochs instead of running to the configured budget of 500. Minibatch ordering is drawn from a seeded `MersenneTwister`, so the presentation order is fixed for a given configuration. Checkpoints are written periodically, which retains both the final model state and the state with the lowest validation loss. The trained checkpoint is distributed with the dataset.

## 7. Evaluation protocol

Training and validation spectra are cut into overlapping windows of 4096 points with a hop of 2048, so adjacent crops share half their width. On a 65536 point spectrum that yields roughly 31 windows. Tiling keeps the input width constant, holds the memory footprint within reach of a single GPU, and exposes the model to local structure from every region of the spectrum. Held out evaluation follows a different protocol. Test spectra are zero filled, apodized, and Fourier transformed in the same way, then evaluated whole at their full 65536 points. Every crop edge in a tiled evaluation introduces a boundary effect, and those effects would otherwise leak into the reported restoration metrics. Keeping the test spectra intact removes that path, at the cost of running the network on much longer inputs, which the padding scheme in section 4 already accommodates. Validation and test data are stratified so that performance is reported separately across the range of conditions in the released datasets instead of as a single aggregate number.
