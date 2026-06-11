# Cross-package validation of `variogram` against GeoStatsFunctions.jl
# (EmpiricalVariogram, Matheron estimator).
#
# Not part of the default CI suite (GeoStatsFunctions is a heavy dependency).
# Run manually in an environment where it is installed, e.g.:
#
#     julia -e 'using Pkg; Pkg.activate(temp=true);
#               Pkg.add(["GeoStatsFunctions","GeoTables","Meshes","Unitful"]);
#               Pkg.develop(path="."); include("test/crosscheck_variogram_geostats.jl")'
#
# Binning is aligned (same maxlag = max_frac × dmax, same number of lags), so
# agreement is at machine precision — both packages compute the same Matheron
# semivariance γ̂(h) = ½·mean[(êᵢ−êⱼ)²] per distance bin.
#
# Verified 2026-06-11: max per-bin rel. diff = 6.0e-16, identical pair counts
# (20/20 bins, n = 400 planar points).

using SpatialHAC, GeoStatsFunctions, GeoTables, Meshes, Unitful, Random

rng = Xoshiro(13)
n = 400
x = 100.0 .* rand(rng, n)            # planar coordinates (km-like unit)
y = 100.0 .* rand(rng, n)
e = randn(rng, n)
per = ones(Int, n)                   # single period

NB = 20
MAXFRAC = 2 / 3

vg = variogram(e, x, y, per; nbins = NB, max_frac = MAXFRAC, distance = :euclidean)
dmax = maximum(hypot(x[i] - x[j], y[i] - y[j]) for i in 1:n for j in (i+1):n)
hmax = MAXFRAC * dmax

gt = georef((z = e,), [(x[i], y[i]) for i in 1:n])
γ = EmpiricalVariogram(gt, :z; nlags = NB, maxlag = hmax, estimator = :matheron)

maxrel = 0.0
nshared = 0
for b in 1:NB
    (vg.n_pairs[b] > 0 && γ.counts[b] > 0) || continue
    global maxrel = max(maxrel, abs(vg.value[b] - γ.ordinates[b]) / abs(γ.ordinates[b]))
    global nshared += 1
end

println("shared non-empty bins: $nshared / $NB")
println("max per-bin rel. diff: $maxrel")
println("pair counts equal: ", vg.n_pairs == Int.(γ.counts))
maxrel < 1e-12 || error("GeoStats crosscheck exceeds 1e-12: $maxrel")
vg.n_pairs == Int.(γ.counts) || error("pair counts differ")
println("GeoStatsFunctions crosscheck PASSED (machine precision)")
