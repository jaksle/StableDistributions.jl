module StableDistributions

using Random
import Random: default_rng, rand!, SamplerRangeInt
using Distributions
using SpecialFunctions
import QuadGK: quadgk
import Distributions: @check_args, @distr_support, @quantile_newton,
    params, shape, location, scale, minimum, maximum,
    convert, convolve, +, *, partype, mean, var, skewness, kurtosis,
    quantile_newton, cquantile_newton, invlogcdf_newton, invlogccdf_newton,
    mgf, cf, pdf, logpdf, cdf, quantile, fit


export Stable,
    rand,
    params, shape, location, scale, minimum, maximum,
    convert, convolve, +, *, partype, mean, var, skewness, kurtosis,
    mgf, cf, pdf, logpdf, cdf, mgf, quantile, fit, fit_quantile

include("stable.jl")
include("conversion.jl")
include("convolution.jl")

end
