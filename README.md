# StableDistributions.jl

Generation and estimation of the class of stable distributions (see on [Wikipedia](https://en.wikipedia.org/wiki/Stable_distribution)) in Julia. It fully complies with the interface of `Distributions` package.

## Defining stable distribution

This package uses the so-called type-1 parametrisation determined by:
- stability index 0 < `α` ≤ 2,
- skewness parameter -1 ≤ `β` ≤ 1,
- scale 0 < `σ`,
- location `μ`

Such distribution is uniquely characterised by its characteristic function (Fourier transform of its pdf)
```math
\varphi(t; \alpha, \beta, \sigma, \mu) = \exp\big(\mathrm i t\mu -|\sigma t|^\alpha(1-\mathrm i\beta\mathrm{sgn}(t)\Phi(t))\big)
```
with $\Phi(t) = -\frac{2}{\pi}\log|t|$ for α = 1 or $\Phi(t) = \tan(\pi\alpha/2)$ for α ≠ 1.

To construct stable distribution one can use:
- `Stable(α) ` for standard symmetric α-Stable distribution equivalent to Stable(α, 0, 1, 0),
- `Stable(α, β)` for standard α-Stable distribution with skewness parameter β equivalent to S(α, β, 1, 0),
- `Stable(α, β, σ, μ)` for any parameters.

## Generating values from stable distribution

One can use standard functions `rand(d::Stable)` for one value or `rand(d::Stable, shape)` for a series of values formatted with a given shape.

## Utility functions

One can multiply stable distributions by scalars and add them, e.g. `2Stable(1.5) + 3` is a valid code which returns `Stable(1.5, 0, 2, 3)`. Function `convolve` is used to convolve two stable distributions with the same α, which is a distribution of the sum of two indepedent stable variables with such distributions. For example `convolve(Stable(0.5,1), Stable(0.5,0,2,1))` gives approximately `Stable(0.5, 0.41, 5.83, 1)`. Function `support` returns support of a given stable distribution, which can be half-bounded for skewed stable distribution with β = 1 or β = -1 and α < 1.

[![Build Status](https://github.com/jaksle/StableDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jaksle/StableDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
