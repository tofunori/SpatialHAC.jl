"""
    SpatialHAC

Spatial HAC (Conley) standard errors for fitted linear mixed models at scale.

Computes the Liang-Zeger GLS sandwich for the fixed effects of a
`MixedModels.LinearMixedModel`, with a distance-decaying (Bartlett) kernel over
same-period spatial pairs:

    Var(β̂) = B · M · B
    B = (X'Ω⁻¹X)⁻¹                       (GLS bread)
    M = A'(K ⊙ êê')A = Σᵢⱼ K(dᵢⱼ)·rᵢ·rⱼ'  (spatial meat)
    A = Ω⁻¹X ;  rᵢ = Aᵢ·êᵢ ;  ê = y − Xβ̂  (marginal residuals)
    Ω = I + WW',  W = ZΛ                  (from the fitted model)

Point estimates and Wald SEs are unchanged; only a spatially robust covariance
of the fixed effects is added. Ω⁻¹ is applied via the Woodbury identity with
one sparse Cholesky of `I + W'W`, so the method scales to millions of rows for
nested random-effect designs.

Pitfall this package exists to avoid: applying Ω⁻¹ to the *residuals*
(`sᵢ = xᵢ·(Ω⁻¹ê)ᵢ`) instead of the *design* is algebraically different and
invalid for a GLS estimator — and the error is undetectable in the OLS limit,
where both constructions coincide. The test suite anchors the formula against
an independent dense definitional sandwich and the textbook OLS Conley.

See `vcov_conley` for the main entry point.
"""
module SpatialHAC

using LinearAlgebra, SparseArrays, Statistics, Random
using MixedModels
using MixedModels: varest
using MixedModels: CoefTable
using SpecialFunctions: erfc, erfcinv
import StatsAPI
using StatsAPI: coef, response, nobs, coefnames

export vcov_conley, ConleyResult, suggest_cutoff, CovariogramResult,
       vcov_cluster, ClusterResult,
       covariogram, variogram, SpatialDiagnostic

const EARTH_RADIUS_KM = 6371.0088

"""
    ConleyResult

Result for one cutoff: fields `cutoff` (km), `vcov` (p×p), `se` (Vector),
`coef` (point estimates), `names` (coefficient names), `n_pairs`,
`min_eig` (smallest pre-floor eigenvalue), `floored::Bool`, `kernel::Symbol`,
`distance::Symbol`.

Implements the `StatsAPI` accessors `vcov`, `stderror`, `coefnames` and
`coeftable`; `show` prints a coefficient table.
"""
struct ConleyResult
    cutoff::Float64
    vcov::Matrix{Float64}
    se::Vector{Float64}
    coef::Vector{Float64}
    names::Vector{String}
    n_pairs::Int
    min_eig::Float64
    floored::Bool
    kernel::Symbol
    distance::Symbol
end

"""
    ClusterResult

Result of `vcov_cluster`: fields `vcov` (p×p), `se` (Vector), `coef` (point
estimates), `names` (coefficient names), `n_clusters`, `type`
(`:CR0`/`:CR1`/`:CR1S`), `dof` (residual degrees of freedom, N−p).

Implements the `StatsAPI` accessors `vcov`, `stderror`, `coefnames` and
`coeftable`; `show` prints a coefficient table.
"""
struct ClusterResult
    vcov::Matrix{Float64}
    se::Vector{Float64}
    coef::Vector{Float64}
    names::Vector{String}
    n_clusters::Int
    type::Symbol
    dof::Int
end

"""
    _kernel_fun(kernel::Symbol) -> Function

Map a kernel name to its weight `K(u)`, `u = d/cutoff ∈ [0, 1)`:
`:bartlett` `1−u` (K₁), `:bartlett2` `(1−u)²` (K₂), `:uniform` `1`,
`:epanechnikov` `1−u²`. All satisfy `K(0)=1`.
"""
function _kernel_fun(kernel::Symbol)
    kernel === :bartlett     && return u -> 1.0 - u
    kernel === :bartlett2    && return u -> (1.0 - u)^2
    kernel === :uniform      && return u -> 1.0
    kernel === :epanechnikov && return u -> 1.0 - u^2
    throw(ArgumentError("unknown kernel :$(kernel) " *
        "(use :bartlett, :bartlett2, :uniform or :epanechnikov)"))
end

# Kernels whose Gram matrix is positive semi-definite for points in ℝ² (Schoenberg
# class P₂; Kelejian & Prucha 2007 via Golubov 1981). The linear Bartlett K₁ is in
# P₁ only — PSD-guaranteed in 1-D, not in 2-D — hence the eigenvalue flooring.
const _PSD2D_KERNELS = (:bartlett2,)

"""Great-circle distance in km (haversine)."""
function haversine_km(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)::Float64
    phi1 = deg2rad(lat1); phi2 = deg2rad(lat2)
    dphi = phi2 - phi1
    dlmb = deg2rad(lon2 - lon1)
    a = sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlmb / 2)^2
    return 2 * EARTH_RADIUS_KM * asin(min(1.0, sqrt(a)))
end

"""Planar Euclidean distance (coordinates in the same unit as the cutoff)."""
euclidean_dist(x1::Float64, y1::Float64, x2::Float64, y2::Float64)::Float64 =
    hypot(x1 - x2, y1 - y2)

"""Map a `distance` name to its function (validates the name)."""
function _distance_fun(distance::Symbol)
    distance === :haversine && return haversine_km
    distance === :euclidean && return euclidean_dist
    throw(ArgumentError("distance must be :haversine or :euclidean"))
end

"""
    scaled_re_matrix(m::LinearMixedModel) -> SparseMatrixCSC

Reconstruct `W = Z·Λ` from the public MixedModels API: `sparse(rt)` is the raw
`Z` block (n × nlev·S) and `rt.λ` the S×S relative covariance factor, so
`Wᵢ = Zᵢ · kron(I_nlev, λᵢ)`.
"""
function scaled_re_matrix(m)::SparseMatrixCSC{Float64,Int}
    blocks = SparseMatrixCSC{Float64,Int}[]
    for rt in m.reterms
        Zi = sparse(rt)
        li = Matrix(rt.λ)
        S = size(li, 1)
        nlev = size(Zi, 2) ÷ S
        push!(blocks, S == 1 ? Zi .* li[1, 1] :
                      Zi * kron(sparse(I, nlev, nlev), sparse(li)))
    end
    return hcat(blocks...)
end

"""
    _gls_bread_scores(m; check_tol=1e-6) -> (bread, S, n, p)

Shared GLS bread `B = (X'Ω⁻¹X)⁻¹` and per-row scores `S[i,:] = (Ω⁻¹X)ᵢ·êᵢ`
for the fixed-effects sandwich, used by both `vcov_conley` and `vcov_cluster`.
Case weights are handled by whitening (`X`, `ê`, `W = ZΛ` scaled by `√w`).
Runtime-validates that the reconstructed `Ω` reproduces `vcov(m)` and errors
otherwise (closing the only silent-corruption channel).
"""
function _gls_bread_scores(m; check_tol::Float64 = 1e-6)
    X = m.X
    n, p = size(X)
    # Case weights: in the weight-whitened space the GLS sandwich is identical
    # once X, ê and W = ZΛ are each scaled by √w. `m.sqrtwts` is √(prior
    # weights), empty for an unweighted fit (→ √w = 1). The self-check below
    # holds only for this scaling, so it also guards against double-scaling.
    sqw = isempty(m.sqrtwts) ? ones(n) : Float64.(m.sqrtwts)
    ehat = sqw .* (response(m) .- X * coef(m))   # √w-scaled marginal residuals
    Xw = sqw .* X                                # √w-scaled design X̃
    W = Diagonal(sqw) * scaled_re_matrix(m)      # √w-scaled W̃ = √w⊙(ZΛ)
    q = size(W, 2)
    F = cholesky(Symmetric(sparse(1.0I, q, q) + W'W))
    OinvX = Xw .- W * (F \ (W' * Xw))            # A = Ω̃⁻¹X̃ (Woodbury, whitened)
    bread = inv(Symmetric(Xw' * OinvX))

    ref = Matrix(vcov(m))
    relerr = maximum(abs.(varest(m) .* bread .- ref)) / maximum(abs.(ref))
    relerr < check_tol || error(
        "SpatialHAC: W reconstruction does not reproduce vcov(m) " *
        "(rel. err = $(relerr)); refusing to compute robust SEs.")

    S = OinvX .* ehat                            # row scores rᵢ = Aᵢ·êᵢ
    return Matrix(bread), S, n, p
end

"""
    vcov_conley(m, lat, lon, period, cutoffs; kernel=:bartlett, check_tol=1e-6)
        -> Vector{ConleyResult}

Spatial-HAC (Conley) covariance of the fixed effects of a fitted
`LinearMixedModel`. Case-weighted fits are supported: in the weight-whitened
space the GLS sandwich is identical once `X`, the residuals and `W = ZΛ` are
scaled by `√w` (done internally); for an unweighted fit `√w = 1`.

Arguments:
- `m`: fitted model.
- `lat`, `lon`: per-row coordinates in degrees, in the model's row order.
- `period`: per-row integer period (e.g. year); pairs are restricted to the
  same period (panel spatial-HAC convention).
- `cutoffs`: one or more cutoff distances, in km for `:haversine` or in the
  coordinate unit for `:euclidean`.
- `kernel`: `:bartlett` (`1−u`, default), `:bartlett2` (`(1−u)²`), `:uniform`
  (`1`) or `:epanechnikov` (`1−u²`), with `u = d/cutoff`.
- `distance`: `:haversine` (default; `lat`/`lon` in degrees, cutoff in km) or
  `:euclidean` (planar `x`/`y` in the same unit as the cutoff, e.g. projected
  metres/km). Use `:euclidean` when coordinates are already projected.

**Positive semi-definiteness.** The estimator is PSD-admissible only when the
kernel's Gram matrix is PSD for the coordinate dimension (Schoenberg class
`Pₚ`; Kelejian & Prucha 2007). For 2-D coordinates this holds for
`:bartlett2` (`(1−u)² ∈ P₂`) but **not** for the linear `:bartlett` (`1−u ∈ P₁`
only), `:uniform` or `:epanechnikov`. When a non-PSD kernel produces a negative
eigenvalue the covariance is eigenvalue-floored (`floored=true`) and a warning
is issued; switch to `:bartlett2` for a guaranteed-PSD estimator.

All cutoffs are computed in a single spatial sweep at the largest cutoff,
threaded over periods. A runtime self-check verifies that the reconstructed
`Ω` reproduces the model's own `vcov(m)` (tolerance `check_tol`) and errors
otherwise, so silently wrong results are impossible.

Caveats: the cutoff truncation makes this a LOWER BOUND on the spatial
correction (report a cutoff-sensitivity curve); the candidate grid is not
antimeridian-aware (pairs straddling ±180° longitude would be missed).
"""
function vcov_conley(m, lat::AbstractVector, lon::AbstractVector,
                     period::AbstractVector, cutoffs::AbstractVector{<:Real};
                     kernel::Symbol = :bartlett, distance::Symbol = :haversine,
                     check_tol::Float64 = 1e-6)
    isempty(cutoffs) && throw(ArgumentError("no cutoffs given"))
    any(c -> c <= 0, cutoffs) && throw(ArgumentError("cutoffs must be > 0"))
    kfun = _kernel_fun(kernel)                    # also validates the kernel name
    distfun = _distance_fun(distance)             # also validates the distance name
    bread, S, n, p = _gls_bread_scores(m; check_tol = check_tol)
    lat = Float64.(lat); lon = Float64.(lon); period = Int.(period)
    length(lat) == n && length(lon) == n && length(period) == n ||
        throw(ArgumentError("coordinate/period vectors must match the model rows (n=$(n))"))
    cf = Float64.(coef(m)); nm = String.(coefnames(m))

    cuts = sort(Float64.(cutoffs))
    nc = length(cuts)
    groups = _period_groups(period)
    group_meats = [[zeros(p, p) for _ in 1:nc] for _ in groups]
    group_counts = [zeros(Int, nc) for _ in groups]
    Threads.@threads for gi in eachindex(groups)
        _accumulate_pairs_multi!(group_meats[gi], group_counts[gi],
                                 S, lat, lon, groups[gi], cuts, kfun, distance, distfun)
    end

    diag_meat = S' * S
    results = ConleyResult[]
    for (ci, cutoff) in enumerate(cuts)
        meat = copy(diag_meat)
        n_pairs = 0
        for gi in eachindex(groups)
            meat .+= group_meats[gi][ci]
            n_pairs += group_counts[gi][ci]
        end
        meat = Symmetric((meat .+ meat') ./ 2)
        V = bread * meat * bread
        V = Symmetric((V .+ V') ./ 2)
        ev = eigen(V)
        min_eig = minimum(ev.values)
        floored = min_eig < 0
        if floored
            V = Symmetric(ev.vectors * Diagonal(max.(ev.values, 0.0)) * ev.vectors')
            if !(kernel in _PSD2D_KERNELS)
                @warn "vcov_conley: kernel :$(kernel) is not PSD-guaranteed for " *
                      "2-D coordinates (it lies in Schoenberg class P₁, not P₂); the " *
                      "covariance was eigenvalue-floored. Use kernel=:bartlett2 for a " *
                      "PSD-guaranteed estimator." maxlog = 1
            end
        end
        Vm = Matrix(V)
        push!(results, ConleyResult(cutoff, Vm, sqrt.(max.(diag(Vm), 0.0)),
                                    cf, nm, n_pairs, min_eig, floored, kernel, distance))
    end
    return results
end

"""
    vcov_cluster(m, cluster_id; type=:CR1, check_tol=1e-6) -> ClusterResult

Cluster-robust covariance of the fixed effects of a fitted `LinearMixedModel`,
on the GLS estimand (same bread `B = (X'Ω⁻¹X)⁻¹` as [`vcov_conley`](@ref)).
The meat sums the per-cluster outer products of the row scores
`rᵢ = (Ω⁻¹X)ᵢ·êᵢ`:

```
M = Σ_g (Σ_{i∈g} rᵢ)(Σ_{i∈g} rᵢ)' ;   Var(β̂) = B · M · B
```

`cluster_id` is a per-row label vector (any type, compared by equality).
Sums are accumulated by cluster — no global N×N matrix is ever formed, so this
scales to millions of rows and thousands of clusters.

`type` sets the small-sample factor (Cameron & Miller 2015; clubSandwich):
- `:CR0`  no correction (Liang & Zeger 1986 sandwich)
- `:CR1`  × `m/(m−1)`               (m = number of clusters)
- `:CR1S` × `(m(N−1))/((m−1)(N−p))` (Stata's default)

Case-weighted fits are supported (scores are formed in the weight-whitened
space). The meat is PSD by construction, so no eigenvalue flooring is needed.
A runtime self-check verifies the reconstructed `Ω` reproduces `vcov(m)`.
"""
function vcov_cluster(m, cluster_id::AbstractVector;
                      type::Symbol = :CR1, check_tol::Float64 = 1e-6)
    type in (:CR0, :CR1, :CR1S) ||
        throw(ArgumentError("type must be :CR0, :CR1 or :CR1S (got :$(type))"))
    bread, S, n, p = _gls_bread_scores(m; check_tol = check_tol)
    length(cluster_id) == n ||
        throw(ArgumentError("cluster_id length must match the model rows (n=$(n))"))

    sums = Dict{Any,Vector{Float64}}()
    @inbounds for i in 1:n
        v = get!(() -> zeros(p), sums, cluster_id[i])
        for k in 1:p
            v[k] += S[i, k]
        end
    end
    G = length(sums)
    meat = zeros(p, p)
    for sg in values(sums)
        meat .+= sg * sg'
    end

    factor = type === :CR0  ? 1.0 :
             type === :CR1  ? G / (G - 1) :
                              (G * (n - 1)) / ((G - 1) * (n - p))   # :CR1S
    V = bread * (factor .* meat) * bread
    V = Symmetric((V .+ V') ./ 2)
    Vm = Matrix(V)
    return ClusterResult(Vm, sqrt.(max.(diag(Vm), 0.0)),
                         Float64.(coef(m)), String.(coefnames(m)), G, type, n - p)
end

# ---- StatsAPI accessors + display -------------------------------------------

const _RobustResult = Union{ConleyResult,ClusterResult}

StatsAPI.vcov(r::_RobustResult) = r.vcov
StatsAPI.stderror(r::_RobustResult) = r.se
StatsAPI.coef(r::_RobustResult) = r.coef
StatsAPI.coefnames(r::_RobustResult) = r.names

"""
    coeftable(r; level=0.95) -> CoefTable

Coefficient table for a robust result: point estimates (unchanged from the
model), robust SEs, Wald `z = est/se`, two-sided normal p-values, and a
`level` confidence interval.
"""
function StatsAPI.coeftable(r::_RobustResult; level::Real = 0.95)
    est = r.coef; se = r.se
    z = est ./ se
    p = erfc.(abs.(z) ./ sqrt(2))                 # two-sided standard-normal
    zc = -sqrt(2) * erfcinv(1 + level)            # e.g. 1.96 at level=0.95
    lo = est .- zc .* se; hi = est .+ zc .* se
    lvl = round(Int, 100 * level)
    CoefTable([est, se, z, p, lo, hi],
              ["Coef.", "Std. Error", "z", "Pr(>|z|)", "Lower $(lvl)%", "Upper $(lvl)%"],
              r.names, 4, 3)
end

function Base.show(io::IO, mime::MIME"text/plain", r::ConleyResult)
    unit = r.distance === :haversine ? " km" : ""
    dist = r.distance === :haversine ? "" : ", distance = :$(r.distance)"
    println(io, "ConleyResult (spatial-HAC, kernel = :$(r.kernel)$(dist), ",
            "cutoff = $(r.cutoff)$(unit), n_pairs = $(r.n_pairs)",
            r.floored ? ", PSD-floored" : "", ")")
    show(io, mime, coeftable(r))
end

function Base.show(io::IO, mime::MIME"text/plain", r::ClusterResult)
    println(io, "ClusterResult ($(r.type), $(r.n_clusters) clusters, dof = $(r.dof))")
    show(io, mime, coeftable(r))
end

"""Group row indices by period value."""
function _period_groups(period::Vector{Int})::Vector{Vector{Int}}
    d = Dict{Int,Vector{Int}}()
    for (i, y) in enumerate(period)
        push!(get!(d, y, Int[]), i)
    end
    return collect(values(d))
end

"""
Accumulate kernel-weighted cross products for all unordered same-period pairs
within each cutoff in `cuts` (sorted), in one sweep at the largest cutoff. The
coordinate grid only proposes candidates; `distfun` decides the distance. Cell
sizes are set so that every within-cutoff pair lands in an adjacent cell: from
degrees-to-km for `:haversine`, or directly in coordinate units for
`:euclidean`.
"""
function _accumulate_pairs_multi!(meats::Vector{Matrix{Float64}},
                                  counts::Vector{Int}, S::Matrix{Float64},
                                  lat::Vector{Float64}, lon::Vector{Float64},
                                  idx::Vector{Int}, cuts::Vector{Float64},
                                  kfun::F, mode::Symbol, distfun::F2) where {F,F2}
    isempty(idx) && return nothing
    p = size(S, 2)
    nc = length(cuts)
    cmax = cuts[end]
    if mode === :haversine
        dlat = cmax / 110.574
        maxabslat = maximum(abs.(@view lat[idx]))
        coslat = max(cosd(min(maxabslat, 89.0)), 1e-3)
        dlon = cmax / (111.320 * coslat)
    else                                    # :euclidean — coords in cutoff units
        dlat = cmax; dlon = cmax
    end

    grid = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in idx
        key = (floor(Int, lat[i] / dlat), floor(Int, lon[i] / dlon))
        push!(get!(grid, key, Int[]), i)
    end

    si = Vector{Float64}(undef, p)
    for i in idx
        ci = floor(Int, lat[i] / dlat); cj = floor(Int, lon[i] / dlon)
        @inbounds for k in 1:p
            si[k] = S[i, k]
        end
        for a in (ci - 1):(ci + 1), b in (cj - 1):(cj + 1)
            cell = get(grid, (a, b), nothing)
            cell === nothing && continue
            for j in cell
                j <= i && continue
                d = distfun(lat[i], lon[i], lat[j], lon[j])
                d >= cmax && continue
                for ck in 1:nc
                    c = cuts[ck]
                    d >= c && continue
                    w = kfun(d / c)
                    counts[ck] += 1
                    meat = meats[ck]
                    @inbounds for r in 1:p, cc in 1:p
                        meat[r, cc] += w * (si[r] * S[j, cc] + S[j, r] * si[cc])
                    end
                end
            end
        end
    end
    return nothing
end

"""
    CovariogramResult

Result of the residual-covariogram cutoff selector (`suggest_cutoff`):

- `cutoff`: selected bandwidth ς̂ in km (`NaN` if the covariogram never
  reached the tolerance band within the binned domain);
- `bins`: distance-bin centers (km);
- `C`: empirical covariogram `Ĉ(h_b)` per bin (`NaN` for empty bins);
- `n_pairs`: number of same-period pairs per bin;
- `n_used`: rows retained after subsampling;
- `crossed`: whether a zero crossing / tolerance hit was found.
"""
struct CovariogramResult
    cutoff::Float64
    bins::Vector{Float64}
    C::Vector{Float64}
    n_pairs::Vector{Int}
    n_used::Int
    crossed::Bool
end

"""
    suggest_cutoff(m::LinearMixedModel, lat, lon, period; kwargs...) -> CovariogramResult
    suggest_cutoff(ehat::AbstractVector, lat, lon, period; kwargs...) -> CovariogramResult

Data-driven Conley cutoff via the covariogram-range selector of Lehner (2026,
arXiv:2603.03997): bin the products `êᵢ·êⱼ` of same-period residual pairs by
great-circle distance, and select the first bin center at which the empirical
covariogram `Ĉ(h)` crosses zero (tolerance `eta`, default 0). The selected
range is used directly as the kernel cutoff — no kernel-dependent rescaling.

For a fitted model the marginal residuals `ê = y − Xβ̂` are used — the same
residuals that enter the `vcov_conley` meat.

Keyword arguments (defaults follow the paper's recommendations):
- `nbins = 150`: number of distance bins (paper: 100-200);
- `max_frac = 2/3`: bins span `max_frac` × the maximum same-period
  inter-point distance;
- `eta = 0.0`: tolerance band; `0` selects the first sign change of `Ĉ`;
- `max_points = 10_000`: rows are subsampled (without replacement, seeded
  `rng`) above this size — the covariogram is O(n²) in pairs. This is our
  adaptation for large panels; the paper benchmarks n ≈ 10⁴ without
  subsampling;
- `rng = Xoshiro(2026)`: RNG for the subsample (fixed seed → reproducible).

Adaptations vs. the paper (which is cross-sectional OLS): pairs are
restricted to the same `period` (the same convention as `vcov_conley`), and
residuals are the mixed-model marginal residuals.

Always inspect the returned covariogram and report a cutoff-sensitivity
curve around ς̂: Lehner shows SE magnitude is inverse-U in the bandwidth, so
neither tiny nor huge cutoffs are conservative.
"""
function suggest_cutoff(m::LinearMixedModel, lat::AbstractVector,
                        lon::AbstractVector, period::AbstractVector; kwargs...)
    ehat = response(m) .- m.X * coef(m)
    return suggest_cutoff(ehat, lat, lon, period; kwargs...)
end

function suggest_cutoff(ehat::AbstractVector, lat::AbstractVector,
                        lon::AbstractVector, period::AbstractVector;
                        nbins::Int = 150, max_frac::Float64 = 2 / 3,
                        eta::Float64 = 0.0, max_points::Int = 10_000,
                        rng::AbstractRNG = Xoshiro(2026),
                        distance::Symbol = :haversine)
    eta >= 0 || throw(ArgumentError("eta must be >= 0"))
    centers, C, counts, n_used = _binned_pairs((a, b) -> a * b, ehat, lat, lon, period;
        nbins = nbins, max_frac = max_frac, max_points = max_points,
        rng = rng, distance = distance)

    # Lehner Eq. (5): ς̂ = min{ h_b : |Ĉ(h_b)| ≤ η }; with η = 0 this is the
    # first sign change of Ĉ relative to its short-range sign.
    cutoff = NaN
    crossed = false
    s0 = 0.0
    for b in 1:nbins
        counts[b] == 0 && continue
        c = C[b]
        if abs(c) <= eta || (s0 != 0.0 && sign(c) != s0) || (s0 == 0.0 && c <= 0)
            cutoff = centers[b]
            crossed = true
            break
        end
        s0 == 0.0 && (s0 = sign(c))
    end
    return CovariogramResult(cutoff, centers, C, counts, n_used, crossed)
end

# ---- shared pair-binning + spatial diagnostics ------------------------------

"""
    _binned_pairs(statfn, ehat, lat, lon, period; kwargs...) -> (h, value, counts, n_used)

Bin a per-pair statistic `statfn(êᵢ, êⱼ)` over same-period pairs by distance.
`statfn = (a,b)->a*b` gives the covariogram; `(a,b)->0.5*(a-b)^2` gives the
semivariogram. Rows are subsampled (seeded) above `max_points` (pairs are
O(n²)). Shared by `suggest_cutoff`, `covariogram` and `variogram`.
"""
function _binned_pairs(statfn::SF, ehat::AbstractVector, lat::AbstractVector,
                       lon::AbstractVector, period::AbstractVector;
                       nbins::Int, max_frac::Float64, max_points::Int,
                       rng::AbstractRNG, distance::Symbol) where {SF}
    n = length(ehat)
    length(lat) == n && length(lon) == n && length(period) == n ||
        throw(ArgumentError("residual/coordinate/period vectors must have equal length"))
    nbins >= 2 || throw(ArgumentError("nbins must be >= 2"))
    0 < max_frac <= 1 || throw(ArgumentError("max_frac must be in (0, 1]"))
    distfun = _distance_fun(distance)

    e = Float64.(ehat); la = Float64.(lat); lo = Float64.(lon); per = Int.(period)
    use = n <= max_points ? collect(1:n) : sort!(randperm(rng, n)[1:max_points])
    groups = _period_groups(per[use])
    for g in groups                       # back to original row indices
        @inbounds for k in eachindex(g)
            g[k] = use[g[k]]
        end
    end

    dmax = 0.0                            # pass 1: max same-period distance
    for g in groups, ii in eachindex(g), jj in (ii + 1):length(g)
        i = g[ii]; j = g[jj]
        d = distfun(la[i], lo[i], la[j], lo[j])
        d > dmax && (dmax = d)
    end
    dmax > 0 || throw(ArgumentError(
        "no same-period pair with positive distance; cannot bin"))

    hmax = max_frac * dmax
    width = hmax / nbins
    sums = zeros(nbins); counts = zeros(Int, nbins)
    for g in groups, ii in eachindex(g), jj in (ii + 1):length(g)
        i = g[ii]; j = g[jj]
        d = distfun(la[i], lo[i], la[j], lo[j])
        d >= hmax && continue
        b = min(nbins, floor(Int, d / width) + 1)
        sums[b] += statfn(e[i], e[j])
        counts[b] += 1
    end
    centers = [(b - 0.5) * width for b in 1:nbins]
    value = [counts[b] > 0 ? sums[b] / counts[b] : NaN for b in 1:nbins]
    return centers, value, counts, length(use)
end

"""
    SpatialDiagnostic

Empirical spatial diagnostic of regression residuals (`covariogram` or
`variogram`): fields `kind` (`:covariogram` or `:semivariogram`), `h` (bin
centers), `value` (`Ĉ(h)` or `γ̂(h)` per bin, `NaN` for empty bins), `n_pairs`
(per bin), `n_used` (rows after subsampling), `distance`.
"""
struct SpatialDiagnostic
    kind::Symbol
    h::Vector{Float64}
    value::Vector{Float64}
    n_pairs::Vector{Int}
    n_used::Int
    distance::Symbol
end

"""
    covariogram(m, lat, lon, period; kwargs...) -> SpatialDiagnostic
    covariogram(ehat, lat, lon, period; kwargs...) -> SpatialDiagnostic

Empirical **covariogram** `Ĉ(h)` of the (marginal) residuals: the binned mean
of products `êᵢ·êⱼ` over same-period pairs at distance `≈ h` (the same quantity
the `suggest_cutoff` selector operates on, exposed as a diagnostic). `Ĉ(h)`
starts near the residual variance at `h→0` and decays toward 0.

Keywords: `nbins=150`, `max_frac=2/3` (bins span `max_frac × max distance`),
`max_points=10_000` (seeded subsample above this; pairs are O(n²)),
`rng=Xoshiro(2026)`, `distance=:haversine` (or `:euclidean`).
"""
covariogram(m::LinearMixedModel, lat, lon, period; kwargs...) =
    covariogram(response(m) .- m.X * coef(m), lat, lon, period; kwargs...)

function covariogram(ehat::AbstractVector, lat::AbstractVector, lon::AbstractVector,
                     period::AbstractVector; nbins::Int = 150, max_frac::Float64 = 2 / 3,
                     max_points::Int = 10_000, rng::AbstractRNG = Xoshiro(2026),
                     distance::Symbol = :haversine)
    h, v, c, nu = _binned_pairs((a, b) -> a * b, ehat, lat, lon, period;
        nbins = nbins, max_frac = max_frac, max_points = max_points,
        rng = rng, distance = distance)
    return SpatialDiagnostic(:covariogram, h, v, c, nu, distance)
end

"""
    variogram(m, lat, lon, period; kwargs...) -> SpatialDiagnostic
    variogram(ehat, lat, lon, period; kwargs...) -> SpatialDiagnostic

Empirical **semivariogram** `γ̂(h) = ½·mean[(êᵢ−êⱼ)²]` of the (marginal)
residuals over same-period pairs at distance `≈ h` (Matheron estimator). Note
this returns `γ` (the semivariogram), not the variogram `2γ` — matching `gstat`
and `GeoStats.jl`. Under second-order stationarity `γ(h) = Ĉ(0) − Ĉ(h)` with
`Ĉ(0)` the residual variance (sill); `γ̂` rises from the nugget toward the sill.

Same keywords as [`covariogram`](@ref).
"""
variogram(m::LinearMixedModel, lat, lon, period; kwargs...) =
    variogram(response(m) .- m.X * coef(m), lat, lon, period; kwargs...)

function variogram(ehat::AbstractVector, lat::AbstractVector, lon::AbstractVector,
                   period::AbstractVector; nbins::Int = 150, max_frac::Float64 = 2 / 3,
                   max_points::Int = 10_000, rng::AbstractRNG = Xoshiro(2026),
                   distance::Symbol = :haversine)
    h, v, c, nu = _binned_pairs((a, b) -> 0.5 * (a - b)^2, ehat, lat, lon, period;
        nbins = nbins, max_frac = max_frac, max_points = max_points,
        rng = rng, distance = distance)
    return SpatialDiagnostic(:semivariogram, h, v, c, nu, distance)
end

function Base.show(io::IO, ::MIME"text/plain", d::SpatialDiagnostic)
    nfull = count(>(0), d.n_pairs)
    println(io, "SpatialDiagnostic ($(d.kind), $(nfull)/$(length(d.h)) non-empty bins, ",
            "n_used = $(d.n_used), distance = :$(d.distance))")
    print(io, "  h ∈ [$(round(minimum(d.h); digits=3)), $(round(maximum(d.h); digits=3))]")
end

end # module
