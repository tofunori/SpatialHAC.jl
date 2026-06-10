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

using LinearAlgebra, SparseArrays, Statistics
using MixedModels
using MixedModels: varest
using StatsAPI: coef, response, nobs

export vcov_conley, ConleyResult

const EARTH_RADIUS_KM = 6371.0088

"""
    ConleyResult

Result for one cutoff: fields `cutoff` (km), `vcov` (p×p), `se` (Vector),
`n_pairs`, `min_eig` (smallest pre-floor eigenvalue), `floored::Bool`.
"""
struct ConleyResult
    cutoff::Float64
    vcov::Matrix{Float64}
    se::Vector{Float64}
    n_pairs::Int
    min_eig::Float64
    floored::Bool
end

"""Great-circle distance in km (haversine)."""
function haversine_km(lat1::Float64, lon1::Float64, lat2::Float64, lon2::Float64)::Float64
    phi1 = deg2rad(lat1); phi2 = deg2rad(lat2)
    dphi = phi2 - phi1
    dlmb = deg2rad(lon2 - lon1)
    a = sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlmb / 2)^2
    return 2 * EARTH_RADIUS_KM * asin(min(1.0, sqrt(a)))
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
    vcov_conley(m, lat, lon, period, cutoffs; check_tol=1e-6) -> Vector{ConleyResult}

Spatial-HAC (Conley/Bartlett) covariance of the fixed effects of a fitted,
unweighted `LinearMixedModel`.

Arguments:
- `m`: fitted model.
- `lat`, `lon`: per-row coordinates in degrees, in the model's row order.
- `period`: per-row integer period (e.g. year); pairs are restricted to the
  same period (panel spatial-HAC convention).
- `cutoffs`: one or more cutoff distances in km (Bartlett kernel `1 − d/c`).

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
                     check_tol::Float64 = 1e-6)
    X = m.X
    n, p = size(X)
    lat = Float64.(lat); lon = Float64.(lon); period = Int.(period)
    length(lat) == n && length(lon) == n && length(period) == n ||
        throw(ArgumentError("coordinate/period vectors must match the model rows (n=$(n))"))
    isempty(cutoffs) && throw(ArgumentError("no cutoffs given"))
    any(c -> c <= 0, cutoffs) && throw(ArgumentError("cutoffs must be > 0 km"))

    # weighted fits change the estimating equations; refuse with a clear error
    # rather than letting the vcov(m) self-check fail confusingly downstream
    if !isempty(m.sqrtwts)
        throw(ArgumentError("vcov_conley supports unweighted models only " *
                            "(this model was fitted with case weights)"))
    end

    ehat = response(m) .- X * coef(m)            # marginal residuals
    W = scaled_re_matrix(m)
    q = size(W, 2)
    F = cholesky(Symmetric(sparse(1.0I, q, q) + W'W))

    OinvX = X .- W * (F \ (W' * X))              # A = Ω⁻¹X (Woodbury)
    bread = inv(Symmetric(X' * OinvX))

    # runtime self-validation against the model's own vcov
    ref = Matrix(vcov(m))
    relerr = maximum(abs.(varest(m) .* bread .- ref)) / maximum(abs.(ref))
    relerr < check_tol || error(
        "vcov_conley: W reconstruction does not reproduce vcov(m) " *
        "(rel. err = $(relerr)); refusing to compute spatial SEs.")

    S = OinvX .* ehat                            # row scores rᵢ = Aᵢ·êᵢ

    cuts = sort(Float64.(cutoffs))
    nc = length(cuts)
    groups = _period_groups(period)
    group_meats = [[zeros(p, p) for _ in 1:nc] for _ in groups]
    group_counts = [zeros(Int, nc) for _ in groups]
    Threads.@threads for gi in eachindex(groups)
        _accumulate_pairs_multi!(group_meats[gi], group_counts[gi],
                                 S, lat, lon, groups[gi], cuts)
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
        end
        Vm = Matrix(V)
        push!(results, ConleyResult(cutoff, Vm, sqrt.(max.(diag(Vm), 0.0)),
                                    n_pairs, min_eig, floored))
    end
    return results
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
Accumulate Bartlett-weighted cross products for all unordered same-period pairs
within each cutoff in `cuts` (sorted), in one sweep at the largest cutoff. The
lat/lon grid only proposes candidates; the haversine distance decides.
"""
function _accumulate_pairs_multi!(meats::Vector{Matrix{Float64}},
                                  counts::Vector{Int}, S::Matrix{Float64},
                                  lat::Vector{Float64}, lon::Vector{Float64},
                                  idx::Vector{Int}, cuts::Vector{Float64})
    isempty(idx) && return nothing
    p = size(S, 2)
    nc = length(cuts)
    cmax = cuts[end]
    dlat = cmax / 110.574
    maxabslat = maximum(abs.(@view lat[idx]))
    coslat = max(cosd(min(maxabslat, 89.0)), 1e-3)
    dlon = cmax / (111.320 * coslat)

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
                d = haversine_km(lat[i], lon[i], lat[j], lon[j])
                d >= cmax && continue
                for ck in 1:nc
                    c = cuts[ck]
                    d >= c && continue
                    w = 1.0 - d / c
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

end # module
