module StableDistributions

using Random
import Random: rand, AbstractRNG
using Distributions
using SpecialFunctions
import QuadGK: quadgk
import Distributions: @check_args, @distr_support,
    params, shape, location, scale, support, minimum, maximum,
    convert, convolve, +, *,
    partype, mean, var, skewness, kurtosis,
    quantile_newton, cquantile_newton, invlogcdf_newton, invlogccdf_newton,
    mgf, cf, pdf, logpdf, cdf, quantile, fit


export Stable,
    rand,
    params, shape, location, scale, support, minimum, maximum,
    convert, convolve, +, *,
    partype, mean, var, skewness, kurtosis,
    mgf, cf, pdf, logpdf, cdf, mgf, fit, fit_quantile,
    quantile, cquantile, invlogcdf, invlogccdf

include("stable.jl")
include("conversion.jl")
include("convolution.jl")

end
