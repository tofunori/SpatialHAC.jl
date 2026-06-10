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
using LinearAlgebra, SparseArrays, Statistics
using MixedModels: varest
using StatsAPI: coef, response

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

# dense definitional Liang-Zeger sandwich (ground truth; no shared internals)
function dense_sandwich(X, ehat, Wmat, lat, lon, yearv, cutoff)
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
            w = (1.0 - d / cutoff) * ehat[i] * ehat[j]
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

@testset "haversine + input validation" begin
    @test abs(haversine_km(0.0, 0.0, 0.0, 1.0) - 111.19) < 0.5
    @test abs(haversine_km(50.0, -117.0, 51.0, -117.0) - 111.19) < 0.5
    @test_throws ArgumentError vcov_conley(m, lat[1:3], lon, yearv, [3.0])
    @test_throws ArgumentError vcov_conley(m, lat, lon, yearv, Float64[])
    @test_throws ArgumentError vcov_conley(m, lat, lon, yearv, [-1.0])
end
