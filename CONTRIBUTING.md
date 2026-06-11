# Contributing to SpatialHAC.jl

Thanks for your interest! Issues and pull requests are welcome.

## Reporting bugs / requesting features

Open a [GitHub issue](https://github.com/tofunori/SpatialHAC.jl/issues) with a
minimal reproducible example (model formula, data shape, the call that fails).
For numerical concerns, please include the `min_eig`/`floored` fields of the
result and, if possible, a comparison against one of the dense reference
implementations in `test/runtests.jl`.

## Pull requests

1. Fork, branch from `main`, and keep the change focused.
2. **Every estimator change must keep the two-anchor validation discipline**
   (see `docs/src/validation.md`): a dense definitional reference test in
   `test/runtests.jl`, and — for new estimators — an external cross-language or
   cross-package check (standalone script in `test/`, like
   `crosscheck_cluster.jl`). A new estimator that cannot be validated against
   an independent external reference will not be merged; this is the package's
   core guarantee.
3. Run the full suite locally: `julia --project -e 'using Pkg; Pkg.test()'`
   (includes Aqua.jl and JET.jl quality gates).
4. Update docstrings, the README validation table, and the docs page that the
   change touches.

## Support

This is a single-maintainer research package developed as part of an MSc
project. Best effort is made to respond to issues within a couple of weeks.
Julia ≥ 1.9, MixedModels 4–5 are supported.
