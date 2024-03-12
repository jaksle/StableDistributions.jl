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
\\varphi(t; \\alpha, \\beta, \\sigma, \\mu) = \\exp\\big(\\mathrm i t\\mu -|\\sigma t|^\\alpha(1-\\mathrm i\\beta\\mathrm{sgn}(t)\\Phi(t))\\big)
```
with ``\\Phi(t) = -\\frac{2}{\\pi}\\log|t|`` for α = 1 or ``\\Phi(t) = \\tan(\\pi\\alpha/2)`` for α ≠ 1.

To construct stable distribution one can use:
- `Stable(α) ` for standard symmetric α-Stable distribution equivalent to Stable(α, 0, 1, 0),
- `Stable(α, β)` for standard α-Stable distribution with skewness parameter β equivalent to S(α, β, 1, 0),
- `Stable(α, β, σ, μ)` for any parameters.

## Generating values stable distribution



[![Build Status](https://github.com/jaksle/StableDistributions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jaksle/StableDistributions.jl/actions/workflows/CI.yml?query=branch%3Amain)
