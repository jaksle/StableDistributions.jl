module StableDistributions

using Random
using Random: rand, AbstractRNG
using Distributions
using SpecialFunctions
using QuadGK: quadgk
using StatsFuns: logexpm1, log1mexp
import Distributions: @check_args, @distr_support,
    params, shape, location, scale, support, minimum, maximum,
    convert, convolve, +, *,
    partype, mean, var, skewness, kurtosis,
    mgf, cf, pdf, logpdf, cdf, fit,
    quantile, cquantile, invlogcdf, invlogccdf


export Stable,
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

end
