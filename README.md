# StableDistributions.jl

[![CI](https://github.com/jaksle/StableDistributions.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jaksle/StableDistributions.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/jaksle/StableDistributions.jl/branch/main/graph/badge.svg?token=E8VHSRSXXR)](https://codecov.io/gh/jaksle/StableDistributions.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Generation and estimation of the class of stable distributions (see on [Wikipedia](https://en.wikipedia.org/wiki/Stable_distribution)) in Julia. It fully complies with the interface of [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) package. The algorithms are based on book John P. Nolan, "Univariate Stable Distributions", Springer 2020 and related publications of John P. Nolan.

## Defining

This package uses the so-called type-1 parametrisation determined by:
- stability index 0 < α ≤ 2,
- skewness parameter -1 ≤ β ≤ 1,
- scale 0 < σ,
- location μ.

Such distribution is uniquely characterised by its characteristic function (Fourier transform of its pdf)
```math
\varphi(t; \alpha, \beta, \sigma, \mu) = \exp\big(\mathrm{i}t\mu -|\sigma t|^\alpha(1-\mathrm i\beta\mathrm{sgn}(t)\Phi(t))\big)
```
with $\Phi(t) = -\frac{2}{\pi}\log|t|$ for α = 1 or $\Phi(t) = \tan(\pi\alpha/2)$ for α ≠ 1.

To construct stable distribution one can use:
- `Stable(α)` for standard symmetric α-stable distribution equivalent to Stable(α, 0, 1, 0),
- `Stable(α, β)` for standard α-stable distribution with skewness parameter β equivalent to Stable(α, β, 1, 0),
- `Stable(α, β, σ, μ)` in general case.

## Sampling

One can use standard functions `rand(d::Stable)` for one value or `rand(d::Stable, shape)` for a series of values formatted with a given shape.

## Features

Probability density function `pdf`, comultative probability function `cdf`, moment generating function `mgf`, characteristic function `cf`, quantiles `quantile`, and few ralated functions are also available. Values of `pdf`, `cdf` and `quantile` are approximate and based on numerical approximations of the corresponding integral representions and additional numerical function inversion for `quantile`. Maximum of the pdf can be numerically approximated using `mode`. For example, `pdf(Stable(1.5), 2)` returns pdf of `Stable(1.5)` at point x = 2 which is approximately `0.084`.

## Transforming

One can multiply stable distributions by scalars and add them, e.g. `2Stable(1.5) + 3` is a valid code which returns `Stable(1.5, 0, 2, 3)`. Function `convolve` is used to convolve two stable distributions with the same α, which is a distribution of the sum of two indepedent stable variables with such distributions. For example `convolve(Stable(0.5,1), Stable(0.5,0,2,1))` gives approximately `Stable(0.5, 0.41, 5.83, 1)`. Function `support` returns support of a given stable distribution, which can be half-bounded for skewed stable distribution with β = 1 or β = -1 and α < 1.

Method `convert` can be used to convert from special cases of Lévy, Cauchy and Gaussian to stable distributions. For example, `convert(Stable,Normal(0,1))` returns equivalent of `Stable(2, 0, √2/2, 1)`.

## Fitting

Given data sample one can find best fitting stable distribution using method `fit(Stable, sample)`. It uses algorithm based on finding parameters which best fit sample characteristic function. This methods is considered to be efficient and quick. Alternatively, one can use `fit_quantile(Stable, sample)` which uses older McCulloch's quantile method in type-0 parametrisation. This method is considered to be generally worse but can be useful for additional checks, e.g. when there are doubts about the reliability of the estimation.
