# SpatialHAC.jl

**Spatial HAC (Conley) standard errors for linear mixed models at scale.**

Computes spatially robust (Conley, 1999) standard errors **directly on the
fixed effects of a fitted `MixedModels.LinearMixedModel`** — the exact
estimand of your mixed model, not an OLS surrogate. Scales to millions of
observations via the Woodbury identity and a single sparse Cholesky.

To our knowledge, no existing package computes spatial-HAC standard errors on
mixed models at this scale: `conleyreg`/`fixest`/`acreg` are OLS/GLM-only, and
`clubSandwich` supports mixed models but is neither spatial nor scalable.

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

res = vcov_conley(m, df.latitude, df.longitude, df.year, [5.0, 15.0, 25.0, 50.0])

for r in res
    @show r.cutoff r.se r.n_pairs r.floored
end
```

Each `ConleyResult` carries `cutoff`, `vcov`, `se`, `n_pairs`, `min_eig`,
`floored`. Point estimates and Wald SEs of the model are unchanged.

All cutoffs are computed in **one spatial sweep** at the largest cutoff
(pair distances computed once), threaded over periods (`julia -t auto`).

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

## Caveats

1. **Lower bound:** cutoff truncation ignores residual correlation beyond the
   cutoff (large-group random effects, synoptic fields) and across periods.
   Report a cutoff-sensitivity curve; cite alongside cluster-robust SEs.
2. The candidate grid is not antimeridian-aware (pairs straddling ±180°
   longitude are missed).
3. Unweighted models only; `Ω̂` treated as known (standard GEE practice).

## References

- Conley, T.G. (1999). GMM estimation with cross sectional dependence.
  *J. Econometrics* 92, 1–45.
- Liang, K.-Y., & Zeger, S.L. (1986). Longitudinal data analysis using
  generalized linear models. *Biometrika* 73, 13–22.
- Kelejian, H.H., & Prucha, I.R. (2007). HAC estimation in a spatial framework.
  *J. Econometrics* 140, 131–154.
- Cameron, A.C., & Miller, D.L. (2015). A practitioner's guide to
  cluster-robust inference. *J. Human Resources* 50, 317–372.

## License

MIT
