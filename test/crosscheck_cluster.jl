# Cross-language validation of `vcov_cluster` against R's clubSandwich.
#
# Not part of the default CI suite (CI runners have no R/clubSandwich). Run
# manually where R + lme4 + clubSandwich are installed:
#
#     julia --project test/crosscheck_cluster.jl
#
# It fits the same LMM in MixedModels.jl and lme4, then compares the
# cluster-robust covariance for CR0/CR1/CR1S. Agreement is limited only by the
# two optimizers' convergence tolerance (θ differs at ~1e-4), so a max relative
# error below 1% confirms the formula and the CR small-sample factors.
#
# Verified 2026-06-11: max rel. err = 1.6e-4 for CR0, CR1, CR1S.

using SpatialHAC, MixedModels, DataFrames, CategoricalArrays, Random, DelimitedFiles
using StatsAPI: coef

function panel(; ng = 15, per = 50, seed = 11)
    rng = Xoshiro(seed); rows = NamedTuple[]
    for g in 1:ng, k in 1:per
        i = (g - 1) * per + k
        x = sin(0.31i); y = 0.4 + 0.3x + 0.5 * sin(g) + 0.6 * randn(rng)
        push!(rows, (g = string("g", g), x = x, y = y))
    end
    DataFrame(rows)
end

df = panel()
csv = tempname() * ".csv"
open(csv, "w") do io
    println(io, "g,x,y")
    for r in eachrow(df); println(io, "$(r.g),$(r.x),$(r.y)"); end
end

dfc = copy(df); dfc.g = categorical(dfc.g)
m = fit(MixedModel, @formula(y ~ 1 + x + (1 | g)), dfc; progress = false)
jdir = mktempdir()
for typ in (:CR0, :CR1, :CR1S)
    writedlm(joinpath(jdir, "julia_$(typ).txt"), vcov_cluster(m, dfc.g; type = typ).vcov)
end

rscript = tempname() * ".R"
open(rscript, "w") do io
    print(io, """
    suppressMessages({library(lme4); library(clubSandwich)})
    df <- read.csv("$csv")
    m <- lmer(y ~ 1 + x + (1|g), data=df, REML=TRUE)
    for (typ in c("CR0","CR1","CR1S")) {
      Vr <- as.matrix(vcovCR(m, cluster=df\$g, type=typ))
      Vj <- as.matrix(read.table(sprintf("$jdir/julia_%s.txt", typ)))
      re <- max(abs(Vr - Vj)) / max(abs(Vr))
      cat(sprintf("%-5s  max relerr = %.3e\\n", typ, re))
      if (re > 1e-2) stop(sprintf("%s crosscheck exceeds 1%%: %.3e", typ, re))
    }
    cat("clubSandwich crosscheck PASSED (<1%)\\n")
    """)
end

if Sys.which("Rscript") === nothing
    @info "Rscript not found — skipping clubSandwich crosscheck"
else
    run(`Rscript $rscript`)
end
