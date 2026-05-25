# Gaussian Markov Random Field package for Julia

This package is based on the book by Håvard Rue and Leonard Held:
[Gaussian Markov Random Fields: Theory and Applications](https://www.routledge.com/Gaussian-Markov-Random-Fields-Theory-and-Applications/Rue-Held/p/book/9781584884323).

## Quick start

Define a first-order intrinsic Gaussian Markov random field (iGMRF) with precision `1.0` on a regular lattice of size `(20, 20)`:

```julia
julia> F = iGMRF(20, 20; order = 1, precision = 1.0)
```

A realization of this field can be generated as follows:

```julia
julia> y = rand(F)
```

The joint log-density can be evaluated as follows:

```julia
julia> logpdf(F, y)
```

The full conditional distributions at each grid cell can be obtained as follows:

```julia
julia> GMRF.full_conditionals(F, y)
```

The log-density of `y[i]` under its full conditional distribution, for each grid cell `i`, can be computed as follows:

```julia
julia> GMRF.full_conditionals_logpdf(F, y)
```

> [!NOTE]
> The value `logpdf(F, y)` is the log-density on the constrained intrinsic subspace, whereas `GMRF.full_conditionals_logpdf(F, y)` returns the local full conditional log-densities of the iGMRF. These conditional densities are useful for Gibbs sampling or pseudo-likelihood calculations, but their product is not equal to the joint density.

The conditional distribution of the unknown grid cells, given that the grid cells in `B` are equal to `x`, can be obtained as follows:

```julia
julia> B = [59; 70; 100; 117; 206; 221; 338; 349; 373; 380]
julia> x = y[B]
julia> GMRF.conditional_distribution(F, B, x)
```
