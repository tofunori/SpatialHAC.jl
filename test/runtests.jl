# SpatialHAC.jl test suite.
#
# The decisive checks anchor the FORMULA, not just the implementation:
#  - production vcov_conley == independent dense definitional sandwich
#    (Omega formed densely, no Woodbury, no grid), at several cutoffs
#  - OLS limit (Omega = I) == textbook Conley computed by brute double loop
# Implementation checks (grid == brute pair count, degenerate cutoff == GLS-HC0,
# haversine, threading) complete the suite.

using Test
using SpatialHAC
using SpatialHAC: haversine_km, scaled_re_matrix
using MixedModels, DataFrames, StatsModels, CategoricalArrays
using LinearAlgebra, SparseArrays, Statistics, Random
using MixedModels: varest
using StatsAPI: coef, response, vcov, stderror, coeftable, coefnames

# ---- synthetic spatial panel -------------------------------------------------
function make_panel(; nper = 380)
    rows = NamedTuple[]
    for yr in (2001, 2002), i in 1:nper
        g  = (i % 28) + 1
        px = (g - 1) * 8 + (i % 8) + 1
        lat = 51.0 + 0.05g + 0.001 * sin(13.7i) + 0.02 * (i % 3)
        lon = -117.0 - 0.07g - 0.0015 * cos(7.3i) - 0.028 * (i % 2)
        x = sin(0.31i + yr)
        y = 0.4 + 0.3x + 0.5 * sin(g) + 0.3 * cos(px) + 0.12 * sin(5i + 2yr)
        push!(rows, (g = string("g", g), pix = string("p", px), yr = yr,
                     lat = lat, lon = lon, x = x, y = y))
    end
    df = DataFrame(rows)
    df.g = categorical(df.g); df.pix = categorical(df.pix)
    return df
end

# dense block cluster-robust sandwich (ground truth; forms Ω densely)
function dense_cluster(X, ehat, Wmat, cl, typ)
    n, p = size(X)
    Omega = Matrix(1.0I, n, n) + Matrix(Wmat * Wmat')
    Ad = Omega \ X
    bread = inv(Symmetric(X' * Ad))
    r = Ad .* ehat                          # rᵢ = (Ω⁻¹X)ᵢ·êᵢ
    sums = Dict{Any,Vector{Float64}}()
    for i in 1:n
        v = get!(() -> zeros(p), sums, cl[i]); v .+= r[i, :]
    end
    G = length(sums); M = zeros(p, p)
    for sg in values(sums); M .+= sg * sg'; end
    f = typ === :CR0 ? 1.0 : typ === :CR1 ? G / (G - 1) :
        (G * (n - 1)) / ((G - 1) * (n - p))
    return bread * (f .* M) * bread
end

# dense definitional Liang-Zeger sandwich (ground truth; no shared internals)
function dense_sandwich(X, ehat, Wmat, lat, lon, yearv, cutoff; kfun = u -> 1.0 - u)
    n, p = size(X)
    Omega = Matrix(1.0I, n, n) + Matrix(Wmat * Wmat')
    Ad = Omega \ X
    bread = inv(Symmetric(X' * Ad))
    Sig = zeros(n, n)
    for i in 1:n
        Sig[i, i] = ehat[i]^2
        for j in (i+1):n
            yearv[i] == yearv[j] || continue
            d = haversine_km(lat[i], lon[i], lat[j], lon[j])
            d >= cutoff && continue
            w = kfun(d / cutoff) * ehat[i] * ehat[j]
            Sig[i, j] = w; Sig[j, i] = w
        end
    end
    return bread * (Ad' * Sig * Ad) * bread
end

df = make_panel()
m = fit(MixedModel, @formula(y ~ 1 + x + (1 | g) + (1 | pix) + (0 + x | g)), df;
        progress = false)
X = m.X; n, p = size(X)
lat = Float64.(df.lat); lon = Float64.(df.lon); yearv = Int.(df.yr)
ehat = response(m) .- X * coef(m)
W = scaled_re_matrix(m)

@testset "W reconstruction reproduces vcov(m)" begin
    q = size(W, 2)
    F = cholesky(Symmetric(sparse(1.0I, q, q) + W'W))
    A = X .- W * (F \ (W' * X))
    myv = varest(m) .* inv(Symmetric(X'A))
    ref = Matrix(vcov(m))
    @test maximum(abs.(myv .- ref)) / maximum(abs.(ref)) < 1e-8
end

@testset "production == dense definitional sandwich (multi-cutoff)" begin
    cuts = [1.5, 3.0, 6.0]
    res = vcov_conley(m, lat, lon, yearv, cuts)
    for (k, c) in enumerate(cuts)
        Vdef = dense_sandwich(X, ehat, W, lat, lon, yearv, c)
        @test maximum(abs.(res[k].vcov .- Vdef)) / maximum(abs.(Vdef)) < 1e-9
    end
end

# Healthy weighted panel: a single variance component on 15 well-populated
# groups with genuine residual noise, so the weighted REML fit stays
# non-degenerate (σ² ≫ 0). A rich RE structure on a near-deterministic panel
# drives σ² → 0 under case weights, which is a fitting pathology, not a
# property of the estimator — so we validate the weighted path on a model the
# data actually support.
function make_weighted_panel(; ng = 15, per = 50, seed = 11)
    rng = Xoshiro(seed); rows = NamedTuple[]
    for yr in (2001, 2002), g in 1:ng, k in 1:per
        i = (g - 1) * per + k
        lat = 51.0 + 0.05g + 0.002 * randn(rng)
        lon = -117.0 - 0.07g - 0.002 * randn(rng)
        x = sin(0.31i + yr)
        y = 0.4 + 0.3x + 0.5 * sin(g) + 0.6 * randn(rng)
        push!(rows, (g = string("g", g), yr = yr, lat = lat, lon = lon, x = x, y = y))
    end
    d = DataFrame(rows); d.g = categorical(d.g); d
end

@testset "weighted model == dense weighted sandwich" begin
    dfw = make_weighted_panel(); nw = nrow(dfw)
    latw = Float64.(dfw.lat); lonw = Float64.(dfw.lon); yw = Int.(dfw.yr)
    rng = Xoshiro(7); wts = rand(rng, nw) .+ 0.5
    mw = fit(MixedModel, @formula(y ~ 1 + x + (1 | g)), dfw;
             wts = wts, progress = false)
    @test varest(mw) > 1e-3                       # guard: fit is non-degenerate

    # independent ground truth: scale the raw design/residuals/W by √w (the
    # weights the user passed) and run the same dense Liang-Zeger sandwich.
    sqw = sqrt.(wts)
    Xw = sqw .* mw.X
    ehatw = sqw .* (response(mw) .- mw.X * coef(mw))
    Ww = Diagonal(sqw) * scaled_re_matrix(mw)

    # the internal self-check (varest·bread ≈ vcov(mw)) must hold for weighted too
    q = size(Ww, 2)
    F = cholesky(Symmetric(sparse(1.0I, q, q) + Ww'Ww))
    A = Xw .- Ww * (F \ (Ww' * Xw))
    myv = varest(mw) .* inv(Symmetric(Xw'A))
    refw = Matrix(vcov(mw))
    @test maximum(abs.(myv .- refw)) / maximum(abs.(refw)) < 1e-8

    cuts = [20.0, 40.0, 80.0]
    res = vcov_conley(mw, latw, lonw, yw, cuts)
    for (k, c) in enumerate(cuts)
        Vdef = dense_sandwich(Xw, ehatw, Ww, latw, lonw, yw, c)
        @test maximum(abs.(res[k].vcov .- Vdef)) / maximum(abs.(Vdef)) < 1e-9
    end
end

@testset "unit weights == unweighted fit" begin
    # wts = 1 (√w = 1) must reproduce the unweighted spatial SEs exactly
    dfw = make_weighted_panel(); nw = nrow(dfw)
    latw = Float64.(dfw.lat); lonw = Float64.(dfw.lon); yw = Int.(dfw.yr)
    m0 = fit(MixedModel, @formula(y ~ 1 + x + (1 | g)), dfw; progress = false)
    m1 = fit(MixedModel, @formula(y ~ 1 + x + (1 | g)), dfw;
             wts = ones(nw), progress = false)
    r0 = vcov_conley(m0, latw, lonw, yw, [40.0])[1]
    r1 = vcov_conley(m1, latw, lonw, yw, [40.0])[1]
    @test maximum(abs.(r0.vcov .- r1.vcov)) / maximum(abs.(r0.vcov)) < 1e-6
end

@testset "vcov_cluster == dense block sandwich (CR0/CR1/CR1S)" begin
    # cluster by the RE grouping (g) and by a coarser 3-region grouping
    region = [string("r", (parse(Int, String(df.g[i])[2:end]) - 1) ÷ 10 + 1)
              for i in 1:n]
    for cl in (df.g, region), typ in (:CR0, :CR1, :CR1S)
        res = vcov_cluster(m, cl; type = typ)
        Vd = dense_cluster(X, ehat, W, cl, typ)
        @test maximum(abs.(res.vcov .- Vd)) / maximum(abs.(Vd)) < 1e-9
        @test res.type === typ
        @test res.dof == n - p
    end
end

@testset "singleton clusters == GLS-HC0" begin
    # each row its own cluster → CR0 meat = Σ rᵢrᵢ' = GLS-HC0; matches the
    # degenerate (tiny-cutoff) Conley sandwich diagonal-meat limit.
    rc = vcov_cluster(m, collect(1:n); type = :CR0)
    rh = vcov_conley(m, lat, lon, yearv, [1e-9])[1]
    @test maximum(abs.(rc.vcov .- rh.vcov)) / maximum(abs.(rh.vcov)) < 1e-7
    @test rc.n_clusters == n
end

@testset "vcov_cluster input validation" begin
    @test_throws ArgumentError vcov_cluster(m, df.g; type = :CR2)
    @test_throws ArgumentError vcov_cluster(m, df.g[1:end-1])
end

@testset "StatsAPI accessors + coeftable + show" begin
    res = vcov_conley(m, lat, lon, yearv, [3.0])[1]
    @test vcov(res) == res.vcov
    @test stderror(res) == res.se
    @test coef(res) == coef(m)                     # estimates unchanged
    @test coefnames(res) == coefnames(m)
    ct = coeftable(res)
    @test ct.cols[1] == coef(m)                    # Coef. column == model estimates
    @test ct.cols[2] == res.se
    @test isapprox(ct.cols[3], coef(m) ./ res.se)  # z = est/se
    txt = sprint(show, MIME("text/plain"), res)
    @test occursin("Coef.", txt) && occursin("cutoff", txt)

    rc = vcov_cluster(m, df.g; type = :CR1)
    @test stderror(rc) == rc.se
    @test coeftable(rc).cols[1] == coef(m)
    @test occursin("CR1", sprint(show, MIME("text/plain"), rc))
end

@testset "kernels == dense definitional sandwich" begin
    kerns = [(:bartlett, u -> 1.0 - u), (:bartlett2, u -> (1.0 - u)^2),
             (:uniform, u -> 1.0), (:epanechnikov, u -> 1.0 - u^2)]
    for (kern, kf) in kerns
        res = vcov_conley(m, lat, lon, yearv, [3.0]; kernel = kern)[1]
        Vdef = dense_sandwich(X, ehat, W, lat, lon, yearv, 3.0; kfun = kf)
        @test maximum(abs.(res.vcov .- Vdef)) / maximum(abs.(Vdef)) < 1e-9
        @test res.kernel === kern
    end
    @test_throws ArgumentError vcov_conley(m, lat, lon, yearv, [3.0]; kernel = :gaussian)
end

@testset "kernel PSD: K₂ guaranteed in 2-D, uniform can fail" begin
    # 3-point chain: uniform Gram = [1 1 0; 1 1 1; 0 1 1] has eigenvalue 1−√2 < 0
    xs = [0.0, 1.0, 2.0]; c = 1.5
    Km(kf) = [(d = abs(xs[i] - xs[j]); d < c ? kf(d / c) : 0.0) for i in 1:3, j in 1:3]
    @test minimum(eigen(Symmetric(Km(u -> 1.0))).values) < -1e-6        # uniform fails
    # K₂ Gram is PSD for an arbitrary 2-D cloud (Schoenberg class P₂)
    rng = Xoshiro(3); P = [(randn(rng), randn(rng)) for _ in 1:60]; cc = 1.5
    K2 = [(d = hypot(P[i][1] - P[j][1], P[i][2] - P[j][2]); d < cc ? (1 - d / cc)^2 : 0.0)
          for i in 1:60, j in 1:60]
    @test minimum(eigen(Symmetric(K2)).values) > -1e-8                  # K₂ PSD
end

@testset "OLS limit == textbook Conley (formula anchor)" begin
    CUT = 3.0
    beta_ols = (X'X) \ (X' * df.y)
    e_ols = df.y .- X * beta_ols
    breadO = inv(Symmetric(X'X))
    meatO = zeros(p, p)
    for i in 1:n
        meatO .+= (X[i, :] * X[i, :]') .* e_ols[i]^2
        for j in (i+1):n
            yearv[i] == yearv[j] || continue
            d = haversine_km(lat[i], lon[i], lat[j], lon[j])
            d >= CUT && continue
            w = (1.0 - d / CUT) * e_ols[i] * e_ols[j]
            meatO .+= w .* (X[i, :] * X[j, :]' .+ X[j, :] * X[i, :]')
        end
    end
    Vbrute = breadO * meatO * breadO
    Wempty = spzeros(n, 1)                      # Omega = I
    Vdef = dense_sandwich(X, e_ols, Wempty, lat, lon, yearv, CUT)
    @test maximum(abs.(Vdef .- Vbrute)) / maximum(abs.(Vbrute)) < 1e-9
end

@testset "degenerate cutoff == GLS-HC0" begin
    res0 = vcov_conley(m, lat, lon, yearv, [1e-9])[1]
    q = size(W, 2)
    F = cholesky(Symmetric(sparse(1.0I, q, q) + W'W))
    A = X .- W * (F \ (W' * X))
    breadG = inv(Symmetric(X'A))
    hc0 = breadG * (A' * (A .* (ehat .^ 2))) * breadG
    @test maximum(abs.(res0.vcov .- hc0)) / maximum(abs.(hc0)) < 1e-9
    @test res0.n_pairs == 0
end

@testset "grid pair count == brute force" begin
    CUT = 3.0
    res = vcov_conley(m, lat, lon, yearv, [CUT])[1]
    nb = 0
    for i in 1:n, j in (i+1):n
        yearv[i] == yearv[j] || continue
        haversine_km(lat[i], lon[i], lat[j], lon[j]) < CUT && (nb += 1)
    end
    @test res.n_pairs == nb
    @test nb > 100
end

# ---- covariogram cutoff selector (Lehner 2026) -------------------------------

# independent brute-force reference: same definition, no shared code path
function brute_covariogram(e, la, lo, per; nbins, max_frac)
    n = length(e)
    dmax = 0.0
    for i in 1:n, j in (i+1):n
        per[i] == per[j] || continue
        dmax = max(dmax, haversine_km(la[i], lo[i], la[j], lo[j]))
    end
    hmax = max_frac * dmax
    width = hmax / nbins
    sums = zeros(nbins); cnt = zeros(Int, nbins)
    for i in 1:n, j in (i+1):n
        per[i] == per[j] || continue
        d = haversine_km(la[i], lo[i], la[j], lo[j])
        d >= hmax && continue
        b = min(nbins, floor(Int, d / width) + 1)
        sums[b] += e[i] * e[j]; cnt[b] += 1
    end
    C = [cnt[b] > 0 ? sums[b] / cnt[b] : NaN for b in 1:nbins]
    centers = [(b - 0.5) * width for b in 1:nbins]
    return centers, C, cnt
end

@testset "covariogram == brute-force reference (no subsample)" begin
    rng = Xoshiro(7)
    ns = 240
    la2 = 51.0 .+ 0.4 .* rand(rng, ns)
    lo2 = -117.0 .- 0.5 .* rand(rng, ns)
    pe2 = [i <= ns ÷ 2 ? 2001 : 2002 for i in 1:ns]
    e2 = randn(rng, ns)
    res = suggest_cutoff(e2, la2, lo2, pe2; nbins = 40, max_points = 10_000)
    centers, Cref, cntref = brute_covariogram(e2, la2, lo2, pe2;
                                              nbins = 40, max_frac = 2 / 3)
    @test res.bins ≈ centers
    @test res.n_pairs == cntref
    for b in 1:40
        cntref[b] == 0 && continue
        @test res.C[b] ≈ Cref[b] atol = 1e-12
    end
    # selection rule replicated independently: first sign change of C
    s0 = 0.0; sel = NaN
    for b in 1:40
        cntref[b] == 0 && continue
        c = Cref[b]
        if (s0 != 0.0 && sign(c) != s0) || (s0 == 0.0 && c <= 0)
            sel = centers[b]; break
        end
        s0 == 0.0 && (s0 = sign(c))
    end
    @test res.cutoff ≈ sel
end

@testset "recovers a known spatial range (spherical GP)" begin
    rng = Xoshiro(11)
    ng = 38                               # 38×38 grid ≈ 1444 points, one period
    la3 = Float64[]; lo3 = Float64[]
    for a in 1:ng, b in 1:ng
        push!(la3, 51.0 + 0.018 * a)      # ~2 km spacing
        push!(lo3, -117.0 - 0.029 * b)
    end
    np = length(la3)
    R = 25.0                              # true range, km
    # spherical covariance (PSD in R^3): exact zero beyond R
    Sig = zeros(np, np)
    for i in 1:np
        Sig[i, i] = 1.0
        for j in (i+1):np
            d = haversine_km(la3[i], lo3[i], la3[j], lo3[j])
            d >= R && continue
            v = 1 - 1.5 * (d / R) + 0.5 * (d / R)^3
            Sig[i, j] = v; Sig[j, i] = v
        end
    end
    L = cholesky(Symmetric(Sig + 1e-8I)).L
    e3 = L * randn(rng, np)
    res = suggest_cutoff(e3, la3, lo3, fill(2001, np);
                         nbins = 60, max_points = np)
    @test res.crossed
    @test 0.5R <= res.cutoff <= 1.6R
    # independent noise → selector picks (near-)first bin, far below R
    res0 = suggest_cutoff(randn(rng, np), la3, lo3, fill(2001, np);
                          nbins = 60, max_points = np)
    @test res0.crossed
    @test res0.cutoff < 0.3R
end

@testset "model method, subsampling, validation" begin
    resm = suggest_cutoff(m, lat, lon, yearv; nbins = 30)
    rese = suggest_cutoff(ehat, lat, lon, yearv; nbins = 30)
    @test resm.cutoff === rese.cutoff && resm.C ≈ rese.C    # marginal residuals
    # subsampling: deterministic (seeded) and uses exactly max_points rows
    r1 = suggest_cutoff(ehat, lat, lon, yearv; nbins = 20, max_points = 200)
    r2 = suggest_cutoff(ehat, lat, lon, yearv; nbins = 20, max_points = 200)
    @test r1.n_used == 200
    @test isapprox(r1.C, r2.C; nans = true)
    @test_throws ArgumentError suggest_cutoff(ehat[1:5], lat, lon, yearv)
    @test_throws ArgumentError suggest_cutoff(ehat, lat, lon, yearv; nbins = 1)
    @test_throws ArgumentError suggest_cutoff(ehat, lat, lon, yearv; eta = -0.1)
    @test_throws ArgumentError suggest_cutoff([1.0, 1.0], [50.0, 50.0],
                                              [-117.0, -117.0], [2001, 2002])
end

@testset "haversine + input validation" begin
    @test abs(haversine_km(0.0, 0.0, 0.0, 1.0) - 111.19) < 0.5
    @test abs(haversine_km(50.0, -117.0, 51.0, -117.0) - 111.19) < 0.5
    @test_throws ArgumentError vcov_conley(m, lat[1:3], lon, yearv, [3.0])
    @test_throws ArgumentError vcov_conley(m, lat, lon, yearv, Float64[])
    @test_throws ArgumentError vcov_conley(m, lat, lon, yearv, [-1.0])
end
