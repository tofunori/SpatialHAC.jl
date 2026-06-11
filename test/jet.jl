# JET static-analysis gate (local/dev — not in the CI suite).
#
# JET tracks Julia's Compiler internals so tightly that it routinely fails to
# precompile on Julia versions it does not yet support, which broke CI on all
# three OS runners. Run this manually on a JET-supported Julia version before
# each release:
#
#     julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("JET");
#               Pkg.develop(path="."); include("test/jet.jl")'
#
# Last verified 2026-06-11: Julia 1.12.5 + JET 0.11.3 — 0 possible errors
# (after the scaled_re_matrix empty-reterms inference fix, commit 8814830).

using JET, SpatialHAC

JET.test_package(SpatialHAC; target_modules = (SpatialHAC,))
println("JET static analysis PASSED")
