# SpatialHAC.jl

*Spatial HAC (Conley) standard errors for linear mixed models at scale.*

SpatialHAC.jl computes spatially robust ([Conley, 1999](https://doi.org/10.1016/S0304-4076(98)00084-0))
standard errors **directly on the fixed effects of a fitted
[`MixedModels.jl`](https://juliastats.org/MixedModels.jl/stable/) model** — the
exact estimand of your mixed model, not an OLS surrogate. It scales to millions
of observations via the Woodbury identity and a single sparse Cholesky
factorization.

To our knowledge, no other package computes spatial-HAC standard errors on
mixed models at this scale: `conleyreg`, `fixest` and `acreg` are OLS/GLM-only,
while `clubSandwich` supports mixed models but is neither spatial nor scalable.

## Why you might need this

When observations close in space share unmodeled influences (weather, smoke
plumes, shared coarse-resolution predictors), model-based (Wald) standard
errors treat near-duplicate observations as independent evidence and become
overconfident — often by large factors. SpatialHAC.jl corrects the uncertainty
**without changing your model or its point estimates**: it replaces only the
covariance of the fixed effects with a kernel-weighted sandwich in which pairs
of nearby same-period observations are allowed to covary.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/tofunori/SpatialHAC.jl")
```

## Quick start

```julia
using MixedModels, SpatialHAC

m = fit(MixedModel,
        @formula(y ~ 1 + x + (1|group) + (1|unit) + (0 + x|group)),
        df)

res = vcov_conley(m, df.latitude, df.longitude, df.year,
                  [5.0, 15.0, 25.0, 50.0])   # cutoffs in km

for r in res
    @show r.cutoff r.se r.n_pairs r.floored
end
```

Each [`ConleyResult`](@ref) carries the cutoff, the robust covariance matrix,
the standard errors, the number of in-range pairs, the smallest pre-floor
eigenvalue and whether eigenvalue flooring was applied. All cutoffs are
computed in **one spatial sweep** at the largest cutoff (each pair distance is
computed once), threaded over periods — start Julia with `julia -t auto`.

## Safety by construction

[`vcov_conley`](@ref) self-validates at runtime: the covariance structure
reconstructed from the fitted model must reproduce the model's own `vcov(m)`
to a 1e-6 relative tolerance, otherwise the function refuses to return
results. This closes the only silent-corruption channel (a wrong `W = ZΛ`
reconstruction, for instance after an internal MixedModels.jl change).

## Validation

| Check | Independent reference | Agreement |
|:--|:--|:--|
| Multi-cutoff production path | Dense definitional sandwich (no Woodbury, no grid) | ~1e-14 |
| OLS limit (`Ω = I`) | Textbook Conley estimator, brute-force double loop | ~1e-15 |
| Degenerate cutoff `→ 0` | GLS-HC0 | ~1e-16 |
| Cross-language | Independent R/`lme4` reimplementation | 1.8e-14 (matched θ) |
| Coverage Monte Carlo | Spatially correlated Gaussian-process errors | ≈90% vs 53% (Wald) |

## Caveats

1. **Lower bound.** Cutoff truncation ignores residual correlation beyond the
   cutoff (large-group random effects, synoptic-scale fields) and across
   periods. Report a cutoff-sensitivity curve and cite these SEs alongside —
   not instead of — cluster-robust SEs.
2. The candidate grid is not antimeridian-aware: pairs straddling ±180°
   longitude are missed (the distances themselves are exact).
3. The estimated covariance structure is treated as known (standard GEE
   practice). Case-weighted fits are supported: the GLS sandwich is computed
   in the weight-whitened space (`X`, residuals and `W = ZΛ` scaled by `√w`).
