# SpatialHAC.jl

[![CI](https://github.com/tofunori/SpatialHAC.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/tofunori/SpatialHAC.jl/actions/workflows/CI.yml)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://tofunori.github.io/SpatialHAC.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Spatial HAC (Conley) standard errors for linear mixed models at scale.**

Computes spatially robust (Conley, 1999) standard errors **directly on the
fixed effects of a fitted `MixedModels.LinearMixedModel`** — the exact
estimand of your mixed model, not an OLS surrogate. Scales to millions of
observations via the Woodbury identity and a single sparse Cholesky.

To our knowledge, no existing package computes spatial-HAC standard errors on
mixed models at this scale: `conleyreg`/`fixest`/`acreg` are OLS/GLM-only, and
`clubSandwich` supports mixed models but is neither spatial nor scalable.

Cluster-robust standard errors (`vcov_cluster`, CR0/CR1/CR1S) are also provided
on the same GLS estimand, for triangulation alongside the spatial correction;
they validate against `clubSandwich` and scale by per-cluster accumulation.
Case-weighted fits are supported throughout.

## The estimator

For the GLS estimator of a fitted LMM (`V = σ²Ω`, `Ω = I + WW'`, `W = ZΛ`):

```
Var(β̂) = B · M · B
B = (X'Ω⁻¹X)⁻¹                        (GLS bread)
M = A'(K ⊙ êê')A = Σᵢⱼ K(dᵢⱼ)·rᵢ·rⱼ'   (kernel-weighted meat)
A = Ω⁻¹X ;  rᵢ = Aᵢ·êᵢ ;  ê = y − Xβ̂   (marginal residuals)
```

Bartlett kernel `K(d) = max(0, 1 − d/cutoff)` (PSD-admissible; Kelejian &
Prucha 2007), pairs restricted to the same period (panel convention),
great-circle distances, negative-eigenvalue flooring if needed.

> **The pitfall this package exists to avoid.** Applying `Ω⁻¹` to the
> *residuals* (`sᵢ = xᵢ·(Ω⁻¹ê)ᵢ`) instead of the *design* is algebraically
> different and invalid for a GLS estimator — and undetectable in the OLS
> limit, where both constructions coincide. The test suite anchors the formula
> against an independent dense definitional sandwich and the textbook OLS
> Conley estimator.

## Usage

```julia
using MixedModels, SpatialHAC

m = fit(MixedModel, @formula(y ~ 1 + x + (1|glacier) + (1|pixel) + (0+x|glacier)), df)

# data-driven cutoff: covariogram-range selector (Lehner 2026, arXiv:2603.03997)
sel = suggest_cutoff(m, df.latitude, df.longitude, df.year)
@show sel.cutoff sel.crossed

res = vcov_conley(m, df.latitude, df.longitude, df.year, [5.0, 15.0, 25.0, 50.0])

for r in res
    @show r.cutoff r.se r.n_pairs r.floored
end

# cluster-robust SEs on the same GLS estimand (triangulation)
cl = vcov_cluster(m, df.glacier; type = :CR1)   # :CR0 / :CR1 / :CR1S
@show cl.se cl.n_clusters cl.type
```

Each `ConleyResult` carries `cutoff`, `vcov`, `se`, `coef`, `names`,
`n_pairs`, `min_eig`, `floored`. Point estimates of the model are unchanged.
Results implement the `StatsAPI` accessors (`vcov`, `stderror`, `coef`,
`coefnames`, `coeftable`) and pretty-print as a coefficient table:

```julia
julia> vcov_cluster(m, df.glacier; type = :CR1)
ClusterResult (CR1, 312 clusters, dof = 18443)
─────────────────────────────────────────────────────────
              Coef.  Std. Error      z  Pr(>|z|)   …
─────────────────────────────────────────────────────────
(Intercept)    …          …         …      …
x              …          …         …      …
─────────────────────────────────────────────────────────
```

All cutoffs are computed in **one spatial sweep** at the largest cutoff
(pair distances computed once), threaded over periods (`julia -t auto`).

`suggest_cutoff` selects the bandwidth as the first zero crossing of the
empirical covariogram of the marginal residuals (same-period pairs, seeded
subsampling above 10k rows). Lehner (2026) shows SE magnitude is **inverse-U**
in the bandwidth — neither tiny nor huge cutoffs are conservative — and that
this selector controls test size where fixed bandwidths fail. Still report a
sensitivity curve around the selected value.

## Safety

`vcov_conley` self-validates at runtime: the reconstructed `Ω` must reproduce
the model's own `vcov(m)` to 1e-6, otherwise it refuses to emit results. This
closes the only silent-corruption channel (a wrong `W = ZΛ` reconstruction,
e.g. after an internal MixedModels change).

## Validation

| Check | Independent reference | Agreement |
|---|---|---|
| Multi-cutoff production path | Dense definitional sandwich (no Woodbury/grid) | ~1e-14 |
| OLS limit (Ω=I) | Textbook Conley, brute double loop | ~1e-15 |
| Degenerate cutoff → 0 | GLS-HC0 | ~1e-16 |
| Cross-language | Independent R/lme4 reimplementation | 1.8e-14 (matched θ) |
| Coverage Monte Carlo | Spatial GP errors, 300 reps | ≈90% vs 53% (Wald) |
| Covariogram selector | Brute-force binned covariogram + selection rule | exact |
| Range recovery | Spherical-GP field, known 25 km range | within [0.5, 1.6]×R |
| Cluster-robust (CR0/CR1/CR1S) | Dense block sandwich | ~1e-15 |
| Cluster-robust cross-language | R `clubSandwich::vcovCR` (lmer) | 1.6e-4 (optimizer-limited) |
| Weighted models | Dense weighted Liang-Zeger sandwich | ~1e-14 |

## Caveats

1. **Lower bound:** cutoff truncation ignores residual correlation beyond the
   cutoff (large-group random effects, synoptic fields) and across periods.
   Report a cutoff-sensitivity curve; cite alongside cluster-robust SEs.
2. The candidate grid is not antimeridian-aware (pairs straddling ±180°
   longitude are missed).
3. `Ω̂` treated as known (standard GEE practice). Case-weighted fits are
   supported (the GLS sandwich is computed in the weight-whitened space).

## References

- Conley, T.G. (1999). GMM estimation with cross sectional dependence.
  *J. Econometrics* 92, 1–45.
- Liang, K.-Y., & Zeger, S.L. (1986). Longitudinal data analysis using
  generalized linear models. *Biometrika* 73, 13–22.
- Kelejian, H.H., & Prucha, I.R. (2007). HAC estimation in a spatial framework.
  *J. Econometrics* 140, 131–154.
- Cameron, A.C., & Miller, D.L. (2015). A practitioner's guide to
  cluster-robust inference. *J. Human Resources* 50, 317–372.
- Lehner, A. (2026). Bandwidth selection for spatial HAC standard errors.
  arXiv:2603.03997.

## License

MIT
