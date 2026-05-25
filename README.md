# Gaussian Markov Random Field package for Julia

[![Build Status Linux and macOS](https://travis-ci.org/jojal5/GMRF.jl.svg?branch=master)](https://travis-ci.org/jojal5/GMRF.jl)

[![codecov.io](http://codecov.io/github/jojal5/GMRF.jl/coverage.svg?branch=master)](http://codecov.io/github/jojal5/GMRF.jl?branch=master)

[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jojal5.github.io/GMRF.jl/stable/)
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://jojal5.github.io/GMRF.jl/latest/)

The package is based on the book of Harvard Rue and Leonard Held: [Gaussian Markov Random Fields--Theory and Applications](https://www.routledge.com/Gaussian-Markov-Random-Fields-Theory-and-Applications/Rue-Held/p/book/9781584884323)

## Quick start

Let define an iGMRF of order 1 with precision 1.0 on a regular lattice of size `(20, 20)`:

```julia
julia> F = iGMRF(20, 20, order = 1, precision = 1)
```

A realization of this field can be obtained as follows:
```julia
julia> y = rand(F)
```

The joint log-density can be obtained as follows:
```julia
julia> logpdf(F,y)
```

The complete conditional distributions for each grid cells can be obtained as follows:
```julia
julia> GMRF.full_conditionals(F,y)
```

The log-density of `y[i]` for each conditional distributions of grid cell `i` can be obtained as follows:
```julia
julia> GMRF.full_conditionals_logpdf(F,y)
```

The conditional distribution of the iGMRF knowing that grid cells in `B` are equal to `x` can be obtained as follows:
```julia
julia> B = [59; 70; 100; 117; 206; 221; 338; 349; 373; 380]
julia> x = y[B]
julia> GMRF.conditional_distribution(F, B, x)
```
