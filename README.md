# A [Julia](http://julialang.org) Interface to [basiclu](https://github.com/ERGO-Code/basiclu)

| **Linux/macOS/Windows/FreeBSD** |
|:-------------------------------:|
| [![build-gh][build-gh-img]][build-gh-url] [![build-cirrus][build-cirrus-img]][build-cirrus-url] |

[build-gh-img]: https://github.com/amontoison/basiclu.jl/workflows/CI/badge.svg?branch=main
[build-gh-url]: https://github.com/amontoison/basiclu.jl/actions
[build-cirrus-img]: https://img.shields.io/cirrus/github/amontoison/basiclu.jl?logo=Cirrus%20CI
[build-cirrus-url]: https://cirrus-ci.com/github/amontoison/basiclu.jl

## How to install

```julia
julia> ]
pkg> add https://github.com/amontoison/basiclu.jl.git
pkg> test basiclu
```

## Content

[basiclu](https://github.com/ERGO-Code/basiclu) implements a sparse LU factorization and an update method that maintains the factorization after column changes to the matrix. It is intended for use in simplex-type algorithms and has been tailored to hypersparse linear programming problems. It provides routines for solving linear systems with a dense or sparse right-hand side. basiclu can be also used to compute a maximum volume basis and the row deficiency of a m-by-n matrix A.
