---
title: 'SpatialHAC.jl: Spatial HAC (Conley) and cluster-robust standard errors for linear mixed models at scale'
tags:
  - Julia
  - spatial statistics
  - mixed models
  - robust standard errors
  - econometrics
  - spatial autocorrelation
authors:
  - name: Thierry Laurent St-Pierre
    orcid: 0009-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Université du Québec à Trois-Rivières, Canada
    index: 1
date: 11 June 2026
bibliography: paper.bib
---

# Summary

Researchers who model spatially indexed observations — satellite pixels,
field plots, weather stations, survey clusters — routinely use linear mixed
models (LMMs) to absorb grouping structure, and routinely face the same
problem afterwards: residual spatial correlation that the random effects do
not capture makes the model's standard errors too small, sometimes by large
factors. The standard remedy in econometrics is the spatial
heteroskedasticity-and-autocorrelation-consistent (HAC) covariance estimator
of @conley:1999, which reweights residual cross-products by a distance kernel.
But existing implementations target ordinary least squares, so mixed-model
practitioners must either refit their model as OLS (changing the estimand) or
go without spatially robust inference.

`SpatialHAC.jl` computes Conley spatial-HAC standard errors **directly on the
fixed effects of a fitted `MixedModels.jl` linear mixed model** [@mixedmodels]
— the generalized-least-squares (GLS) estimand of the model the researcher
actually fit. It scales to millions of observations through the Woodbury
identity and a single sparse Cholesky factorization, supports case-weighted
fits, and complements the spatial estimator with cluster-robust standard
errors (CR0/CR1/CR1S) on the same estimand, a data-driven bandwidth selector
based on the residual covariogram [@lehner:2026], and empirical covariogram
and semivariogram diagnostics.

# Statement of need

Spatial-HAC corrections exist for OLS in several ecosystems — `conleyreg` and
`fixest` in R, `acreg` in Stata — but none of them operates on a mixed model.
Conversely, `clubSandwich` [@clubsandwich] provides cluster-robust (not
spatial) covariance estimators for `lme4` models, with dense per-cluster
adjustment matrices that do not scale to millions of rows. To our knowledge,
no package in R, Python, Stata, or Julia computes spatial-HAC standard errors
on the fixed effects of a fitted mixed model at scale; researchers who need
both random effects and spatial robustness are left without tooling.

The estimator is the GLS sandwich of @liang-zeger:1986 with a
kernel-weighted spatial meat: for a fitted LMM with marginal covariance
$\sigma^2\Omega$, $\Omega = I + WW'$, $W = Z\Lambda$,

$$\widehat{\mathrm{Var}}(\hat\beta) = B\,M\,B, \qquad
B = (X'\Omega^{-1}X)^{-1}, \qquad
M = \sum_{i,j} K(d_{ij})\, r_i r_j', \qquad
r_i = (\Omega^{-1}X)_i\,\hat e_i,$$

with $\hat e = y - X\hat\beta$ the marginal residuals and $K$ a distance
kernel with bandwidth (cutoff) $c$. $\Omega^{-1}$ is applied through the
Woodbury identity, requiring one sparse Cholesky of $I + W'W$; pair
accumulation runs on a spatial grid, threaded over time periods, and all
cutoffs of a sensitivity curve are computed in a single sweep.

A subtle algebraic trap motivates both the package and its validation
discipline: applying $\Omega^{-1}$ to the *residuals* instead of the *design*
($s_i = x_i(\Omega^{-1}\hat e)_i$) yields a plausible-looking but invalid
estimator — and the two constructions coincide exactly in the OLS limit, so
OLS-reduction tests cannot detect the error. The package therefore anchors
every estimator two independent ways (a dense definitional reference
implemented without any of the production shortcuts, and an external
cross-language or cross-package reference), and self-validates at runtime: the
reconstructed $\Omega$ must reproduce the model's own `vcov` or the call
errors rather than returning numbers.

# State of the field

Kernel admissibility follows @kelejian-prucha:2007: a spatial-HAC estimator is
positive semi-definite only for kernels in Schoenberg's class $P_p$ for the
coordinate dimension $p$; in two dimensions this holds for the squared
triangular kernel $(1-u)^2_+$ but not for the linear Bartlett kernel, which is
why the package offers both (plus uniform and Epanechnikov), floors negative
eigenvalues, and warns when flooring triggers on a non-$P_2$ kernel. Bandwidth
selection follows the covariogram-range selector of @lehner:2026, who shows
the bandwidth–SE relationship is inverse-U-shaped, so neither tiny nor huge
cutoffs are conservative. Cluster-robust small-sample corrections follow the
CR taxonomy of @cameron-miller:2015 and the mixed-model application of
@huang:2022. Multi-cutoff reporting implements the sensitivity-curve practice
the spatial-inference literature recommends, and the documentation frames the
truncated estimator explicitly as a lower bound on the correction, to be
reported alongside cluster-robust standard errors.

# Software design

The package is deliberately a *post-estimation* library: it consumes a fitted
`LinearMixedModel` and never refits, which keeps `MixedModels.jl` — a mature,
independently validated optimizer — as the source of $\hat\beta$, $\Omega$,
and the runtime-validation reference. All estimators share one GLS
bread/scores kernel; case weights are handled by whitening ($X$, $\hat e$, $W$
scaled by $\sqrt{w}$), which the runtime self-check guards against
double-scaling. Results implement the `StatsAPI` accessors (`vcov`,
`stderror`, `coeftable`) and print as coefficient tables. A space×time HAC
extension was considered and deliberately **cut**: no practicable external
reference exists for the mixed-model GLS estimand, and the package does not
ship estimators it cannot validate two independent ways.

# Research impact

The package grew out of, and is used in, an MSc research pipeline analyzing
MODIS glacier albedo and wildfire-aerosol deposition across ~360,000 pixels
and 22 years in Western North America (7.9 million pixel-year rows), where
spatial-HAC corrections on mixed-model fire-effect coefficients inflate
standard errors by factors of 2.8–11.8 across 5–50 km cutoffs while preserving
the sign and significance of the deposition effect. The cross-language
validation suite (R/`lme4` re-implementation agreeing to $1.8\times10^{-14}$
at matched parameters; `clubSandwich` agreement to optimizer precision;
`GeoStatsFunctions.jl` variogram agreement to machine precision) doubles as
reproducible reference material for the estimator itself.

# AI usage disclosure

This package was developed with substantial assistance from Anthropic's Claude
(Claude Code with Opus-class models, 2026): code generation and refactoring,
test scaffolding, documentation drafting, and copy-editing of this manuscript.
The author framed the problem, made the core design decisions (post-estimation
architecture, two-anchor validation discipline, the decision to cut the
space×time estimator, kernel and selector choices), reviewed, edited and
validated all AI-assisted outputs, and verified every validation result
reported here against the test suite and external references.

# Acknowledgements

This work was supported by the author's MSc research at the Université du
Québec à Trois-Rivières under the supervision of Christophe Kinnard.
Computations were performed in part on the Narval cluster of the Digital
Research Alliance of Canada.

# References
