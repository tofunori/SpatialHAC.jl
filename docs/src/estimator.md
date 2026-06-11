# The estimator

## Setup

A fitted linear mixed model has marginal covariance ``V = \sigma^2 \Omega`` with

```math
\Omega = I + W W', \qquad W = Z \Lambda,
```

where ``Z`` is the random-effects design and ``\Lambda`` the (relative)
covariance factor implied by the fitted variance components. The fixed effects
solve the generalized-least-squares (GLS) estimating equation

```math
X' \Omega^{-1} (y - X\hat\beta) = 0 .
```

## The spatial sandwich

For the GLS estimator under an arbitrary true covariance
``\operatorname{Cov}(y) = \Sigma``,

```math
\operatorname{Var}(\hat\beta)
  = B \, \big( X'\Omega^{-1} \Sigma\, \Omega^{-1} X \big) \, B,
\qquad B = (X'\Omega^{-1}X)^{-1}.
```

Following [Conley (1999)](https://doi.org/10.1016/S0304-4076(98)00084-0) and
the GEE logic of [Liang & Zeger (1986)](https://doi.org/10.1093/biomet/73.1.13),
``\Sigma`` is estimated by a kernel-weighted outer product of the **marginal**
residuals ``\hat e = y - X\hat\beta``:

```math
\hat\Sigma_{ij} = K(d_{ij})\, \hat e_i \hat e_j ,
```

with a kernel taper ``K(u)``, ``u = d/c`` (default Bartlett ``1-u``; also
``(1-u)^2``, uniform, Epanechnikov), pairs restricted to the same period (the
panel spatial-HAC convention), and great-circle distances. The estimator is
positive-semidefinite only for a kernel whose Gram matrix is PSD in the
coordinate dimension (Schoenberg class ``P_p``;
[Kelejian & Prucha, 2007](https://doi.org/10.1016/j.jeconom.2006.09.005), via
Golubov 1981): in 2-D this holds for ``(1-u)^2 \in P_2`` but not for the linear
Bartlett ``1-u \in P_1`` only — hence the eigenvalue flooring. Writing
``A = \Omega^{-1} X`` and row scores ``r_i = A_i\, \hat e_i``,

```math
M \;=\; A'\,(K \odot \hat e\hat e')\,A \;=\; \sum_{i,j} K(d_{ij})\; r_i r_j',
\qquad
\widehat{\operatorname{Var}}(\hat\beta) = B\, M\, B .
```

``\sigma^2`` cancels between bread and meat. Two limits anchor the construction:
with ``\Omega = I`` it reduces to the textbook OLS Conley estimator, and as
``c \to 0`` it reduces to the GLS analogue of HC0.

## The pitfall this package exists to avoid

It is tempting to form per-observation scores by applying ``\Omega^{-1}`` to
the *residuals*: ``s_i = x_i \,(\Omega^{-1}\hat e)_i``. This yields

```math
X' \operatorname{diag}(\Omega^{-1}\hat e)\; K\; \operatorname{diag}(\Omega^{-1}\hat e)\, X,
```

which is **not** ``A'(K \odot \hat e\hat e')A`` and is invalid for a GLS
estimator — in our application it inflated standard errors by an absurd ×60 and
flipped a strongly significant coefficient to non-significance. Crucially, the
two constructions **coincide in the OLS limit** (``\Omega = I``), so OLS-based
reduction tests cannot detect the error. The test suite therefore anchors the
formula against an independent *dense definitional* sandwich (``\Omega`` formed
explicitly, no Woodbury, no spatial grid) and against an independent R/`lme4`
reimplementation.

## Scalability

``\Omega^{-1}`` is never formed. With ``M_q = I + W'W`` (size ``q \times q``,
``q`` = total number of random effects),

```math
\Omega^{-1} u = u - W\, M_q^{-1} W' u
```

(Woodbury), requiring a single sparse Cholesky of ``M_q``. For nested
random-effect designs (e.g. pixel ⊂ glacier ⊂ region) the factor has nearly no
fill-in: on a real 7.8M-row panel with ``q \approx 425{,}000``, the
factorization takes under a second and ~33 MB. Pair search uses an exact
spatial-grid candidate filter (the haversine distance decides membership), all
cutoffs are accumulated in a single sweep at the largest cutoff, and periods
are processed in parallel threads.

## Data-driven cutoff selection

[`suggest_cutoff`](@ref) implements the covariogram-range selector of Lehner
(2026, arXiv:2603.03997): bin the products ``\hat e_i \hat e_j`` of
same-period residual pairs by distance (default ~150 bins spanning two-thirds
of the maximum inter-point distance), and select the first bin center at which
the empirical covariogram

```math
\hat C(h) = \frac{1}{|N(h)|} \sum_{(i,j) \in N(h)} \hat e_i \hat e_j
```

crosses zero (Lehner's Eq. 5 with tolerance ``\eta = 0``). The selected range
is used **directly** as the kernel cutoff — no kernel-dependent rescaling. In
Lehner's Monte Carlo this selector controls the false-positive rate at or near
the nominal level and outperforms every fixed-bandwidth alternative, with best
size control under Bartlett and Epanechnikov kernels.

Two adaptations relative to the paper (which is cross-sectional OLS): pairs
are restricted to the same period — the same panel convention as
[`vcov_conley`](@ref) — and the residuals are the mixed-model *marginal*
residuals ``\hat e = y - X\hat\beta``, i.e. exactly the residuals that enter
the spatial meat. Above `max_points` rows (default 10,000) the selector works
on a seeded random subsample, since the covariogram is quadratic in pairs; the
paper benchmarks ``n \approx 10^4`` without subsampling.

The selector is validated in the test suite against an independent brute-force
implementation of the binned covariogram and selection rule, and shown to
recover a known spatial range on simulated spherical-covariance Gaussian
fields.

## Practical guidance

- **Choose cutoffs from the data** with [`suggest_cutoff`](@ref), and
  **report a sensitivity curve** over several cutoffs around the selected
  value. Lehner (2026) shows the relationship between bandwidth and SE
  magnitude is *inverse-U shaped*: both too-narrow and too-wide bandwidths
  understate standard errors, so a deliberately huge cutoff is **not** a
  conservative choice.
- **Interpret as a lower bound**: correlation beyond the cutoff and across
  periods is not captured. Present these SEs alongside cluster-robust SEs at a
  level matched to your predictor's assignment scale (Moulton logic).
- The result reports `min_eig`, `floored` and `kernel`. With 2-D coordinates,
  flooring can occur for the linear `:bartlett`, `:uniform` and `:epanechnikov`
  kernels (not PSD-guaranteed in 2-D); use `:bartlett2` for a PSD-guaranteed
  estimator.
