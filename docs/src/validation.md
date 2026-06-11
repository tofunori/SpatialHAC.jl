# Validation

The central claim of this package is not just that the estimators are
*implemented*, but that they are *validated* — every estimator is anchored two
independent ways: against a **dense definitional reference** (the formula
re-implemented from scratch with none of the production shortcuts: no Woodbury
identity, no spatial grid, no per-cluster accumulation), and where possible
against an **external cross-language/cross-package reference** computed by code
we did not write.

## Why this discipline exists

The first (pre-release) version of this estimator applied `Ω⁻¹` to the
*residuals* instead of the *design matrix*. Both constructions coincide exactly
in the OLS limit (`Ω = I`), so OLS-reduction tests and internal-consistency
tests (grid == brute force) all passed — while mixed-model standard errors were
inflated by a factor of ~60. Only checks that anchor the **formula** — a dense
definitional sandwich, and an independent re-implementation in another language
— catch this class of error. Every estimator added since must pass both kinds
of anchor before it ships; anything we cannot validate externally is cut (see
the space×time HAC caveat in the README).

## Validation matrix

Each row is enforced by the test suite (`test/runtests.jl`) on every CI run,
except the cross-language rows, which are standalone scripts (CI runners lack
R / heavy geostatistics dependencies) re-run manually before each release.

| Check | Independent reference | Agreement | Where |
|---|---|---|---|
| Multi-cutoff production path | Dense definitional sandwich (no Woodbury/grid) | ~1e-14 | `runtests.jl` |
| OLS limit (Ω=I) | Textbook Conley, brute double loop | ~1e-15 | `runtests.jl` |
| Degenerate cutoff → 0 | GLS-HC0 | ~1e-16 | `runtests.jl` |
| Cross-language (Conley) | Independent R/lme4 re-implementation | 1.8e-14 (matched θ) | `tests/julia` (engine repo) |
| Coverage Monte Carlo | Spatial GP errors, 300 reps | ≈90% vs 53% (Wald) | engine repo |
| Weighted models | Dense weighted Liang-Zeger sandwich | ~1e-14 | `runtests.jl` |
| Unit weights | Unweighted fit (exact reduction) | <1e-6 | `runtests.jl` |
| Cluster-robust CR0/CR1/CR1S | Dense block sandwich | ~1e-15 | `runtests.jl` |
| Cluster-robust cross-language | R `clubSandwich::vcovCR` (lmer) | 1.6e-4 (optimizer-limited) | `test/crosscheck_cluster.jl` |
| Singleton clusters | GLS-HC0 (exact limit) | <1e-7 | `runtests.jl` |
| Kernels (4) | Dense definitional sandwich per kernel | ~1e-15 | `runtests.jl` |
| K₂ PSD in 2-D / uniform not | Kernel Gram eigenvalues (Schoenberg `P₂`) | exact | `runtests.jl` |
| Euclidean distance | Dense planar sandwich | ~1e-15 | `runtests.jl` |
| Covariogram selector | Brute-force binned covariogram + selection rule | exact | `runtests.jl` |
| Range recovery | Spherical-GP field, known 25 km range | within [0.5, 1.6]×R | `runtests.jl` |
| Covariogram / variogram | Independent brute-force binning | ~1e-12 | `runtests.jl` |
| Semivariogram identity | `γ = m₂ − Ĉ` per bin (exact algebra) | ~1e-12 | `runtests.jl` |
| Variogram cross-package | `GeoStatsFunctions.EmpiricalVariogram` | 6e-16 | `test/crosscheck_variogram_geostats.jl` |
| Software quality | `Aqua.test_all` (11 checks) | pass | `runtests.jl` |
| Static type analysis | `JET.test_package` (local gate; JET is too Julia-version-coupled for CI) | 0 possible errors | `test/jet.jl` |

## The runtime self-check

Beyond the test suite, every call to `vcov_conley` / `vcov_cluster` validates
itself at runtime: the `Ω = I + WW'` reconstructed from the model's random
effects must reproduce the model's own `vcov(m)` to `check_tol` (default 1e-6),

```math
\hat\sigma^2 (X'\Omega^{-1}X)^{-1} \approx \operatorname{vcov}(m),
```

otherwise the call **errors instead of returning numbers**. This closes the
only silent-corruption channel — a wrong `W = ZΛ` reconstruction (e.g. after an
internal MixedModels.jl change) — and, because the identity holds only for the
correct √w-whitening, it also guards the weighted path against double-scaling.

## A worked end-to-end example

```@example validation
using MixedModels, SpatialHAC, DataFrames, CategoricalArrays, Random
using StatsAPI: coeftable

# small synthetic spatial panel
rng = Xoshiro(1)
n = 600
g = rand(rng, 1:15, n)
df = DataFrame(g = categorical(string.("g", g)),
               lat = 51.0 .+ 0.05 .* g .+ 0.01 .* randn(rng, n),
               lon = -117.0 .- 0.07 .* g .+ 0.01 .* randn(rng, n),
               yr = rand(rng, 2001:2003, n),
               x = randn(rng, n))
df.y = 0.4 .+ 0.3 .* df.x .+ 0.5 .* sin.(g) .+ 0.6 .* randn(rng, n)

m = fit(MixedModel, @formula(y ~ 1 + x + (1 | g)), df; progress = false)

# diagnose the residual spatial structure, then select a cutoff
vg  = variogram(m, df.lat, df.lon, df.yr; nbins = 30)
sel = suggest_cutoff(m, df.lat, df.lon, df.yr; nbins = 30)

# spatial-HAC and cluster-robust SEs on the same GLS estimand
conley  = vcov_conley(m, df.lat, df.lon, df.yr, [5.0, 15.0, 30.0])
cluster = vcov_cluster(m, df.g; type = :CR1)

coeftable(cluster)
```
