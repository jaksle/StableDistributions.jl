module StableDistributions

using Random
import Random: rand, AbstractRNG
using Distributions
using SpecialFunctions
using QuadGK: quadgk
using StatsFuns: logexpm1, log1mexp
using PDMats, FillArrays, LinearAlgebra
import Base: size, length
import Distributions: @check_args, @distr_support,
    params, shape, location, scale, support, minimum, maximum,
    convert, convolve, +, *,
    partype, mean, var, skewness, kurtosis,
    mgf, cf, pdf, logpdf, cdf, fit,
    mode, quantile, cquantile, invlogcdf, invlogccdf


export Stable, EllipticStable,
    rand,
    params, shape, location, scale, support, minimum, maximum,
    convert, convolve, +, *,
    partype, mean, var, skewness, kurtosis,
    mgf, cf, pdf, logpdf, cdf, mgf, fit, fit_quantile,
    quantile, cquantile, invlogcdf, invlogccdf


include("quantlinealgs.jl")
include("stable.jl")
include("conversion.jl")
include("convolution.jl")

include("elliptic.jl")

end
