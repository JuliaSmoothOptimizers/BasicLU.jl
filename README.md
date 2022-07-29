# A [Julia](http://julialang.org) Interface to [BASICLU](https://github.com/ERGO-Code/basiclu)

| **Documentation** | **Continuous Integration** | **Coverage** | **DOI** |
|:-----------------:|:--------------------------:|:------------:|:-------:|
| [![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] | [![build-gh][build-gh-img]][build-gh-url] | [![codecov][codecov-img]][codecov-url] | [![DOI](https://zenodo.org/badge/386473966.svg)](https://zenodo.org/badge/latestdoi/386473966)

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaSmoothOptimizers.github.io/BasicLU.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://JuliaSmoothOptimizers.github.io/BasicLU.jl/dev
[build-gh-img]: https://github.com/JuliaSmoothOptimizers/BasicLU.jl/workflows/CI/badge.svg?branch=main
[build-gh-url]: https://github.com/JuliaSmoothOptimizers/BasicLU.jl/actions
[codecov-img]: https://codecov.io/gh/JuliaSmoothOptimizers/BasicLU.jl/branch/main/graph/badge.svg
[codecov-url]: https://app.codecov.io/gh/JuliaSmoothOptimizers/BasicLU.jl

## How to cite

If you use BasicLU.jl in your work, please cite using the format given in [CITATION.bib](CITATION.bib).

## How to install

```julia
julia> ]
pkg> add BasicLU
pkg> test BasicLU
```

## Content

[BASICLU](https://github.com/ERGO-Code/basiclu) implements a sparse LU factorization and an update method that maintains the factorization after column changes to the matrix. It is intended for use in simplex-type algorithms and has been tailored to hypersparse linear programming problems. It provides routines for solving linear systems with a dense or sparse right-hand side. BASICLU can be also used to compute a maximum volume basis and the row deficiency of a m-by-n matrix A.
