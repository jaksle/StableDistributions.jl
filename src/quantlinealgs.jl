

# variants of algorithms from Distributions.jl with added iteration number limits
function quantile_newton(d::ContinuousUnivariateDistribution, p::Real, xs::Real=mode(d), tol::Real=1e-12, lim::Int = 10^6)
    x = xs + (p - cdf(d, xs)) / pdf(d, xs)
    T = typeof(x)
    if 0 < p < 1
        x0 = T(xs)
        n = 1
        while abs(x-x0) > max(abs(x),abs(x0)) * tol && n < lim
            x0 = x
            x = x0 + (p - cdf(d, x0)) / pdf(d, x0)
            n += 1
        end
        return x
    elseif p == 0
        return T(minimum(d))
    elseif p == 1
        return T(maximum(d))
    else
        return T(NaN)
    end
end

function cquantile_newton(d::ContinuousUnivariateDistribution, p::Real, xs::Real=mode(d), tol::Real=1e-12, lim::Int = 10^6)
    x = xs + (ccdf(d, xs)-p) / pdf(d, xs)
    T = typeof(x)
    if 0 < p < 1
        x0 = T(xs)
        n = 1
        while abs(x-x0) > max(abs(x),abs(x0)) * tol  && n < lim
            x0 = x
            x = x0 + (ccdf(d, x0)-p) / pdf(d, x0)
            n += 1
        end
        return x
    elseif p == 1
        return T(minimum(d))
    elseif p == 0
        return T(maximum(d))
    else
        return T(NaN)
    end
end

function invlogcdf_newton(d::ContinuousUnivariateDistribution, lp::Real, xs::Real=mode(d), tol::Real=1e-12, lim::Int = 10^6)
    T = typeof(lp - logpdf(d,xs))
    if -Inf < lp < 0
        x0 = T(xs)
        if lp < logcdf(d,x0)
            x = x0 - exp(lp - logpdf(d,x0) + logexpm1(max(logcdf(d,x0)-lp,0)))
            n = 1
            while abs(x-x0) >= max(abs(x),abs(x0)) * tol  && n < lim
                x0 = x
                x = x0 - exp(lp - logpdf(d,x0) + logexpm1(max(logcdf(d,x0)-lp,0)))
                n += 1
            end
        else
            x = x0 + exp(lp - logpdf(d,x0) + log1mexp(min(logcdf(d,x0)-lp,0)))
            n = 1
            while abs(x-x0) >= max(abs(x),abs(x0))*tol   && n < lim
                x0 = x
                x = x0 + exp(lp - logpdf(d,x0) + log1mexp(min(logcdf(d,x0)-lp,0)))
                n += 1
            end
        end
        return x
    elseif lp == -Inf
        return T(minimum(d))
    elseif lp == 0
        return T(maximum(d))
    else
        return T(NaN)
    end
end

function invlogccdf_newton(d::ContinuousUnivariateDistribution, lp::Real, xs::Real=mode(d), tol::Real=1e-12, lim::Int = 10^6)
    T = typeof(lp - logpdf(d,xs))
    if -Inf < lp < 0
        x0 = T(xs)
        if lp < logccdf(d,x0)
            x = x0 + exp(lp - logpdf(d,x0) + logexpm1(max(logccdf(d,x0)-lp,0)))
            n = 1
            while abs(x-x0) >= max(abs(x),abs(x0)) * tol   && n < lim
                x0 = x
                x = x0 + exp(lp - logpdf(d,x0) + logexpm1(max(logccdf(d,x0)-lp,0)))
                n += 1
            end
        else
            x = x0 - exp(lp - logpdf(d,x0) + log1mexp(min(logccdf(d,x0)-lp,0)))
            n = 1
            while abs(x-x0) >= max(abs(x),abs(x0)) * tol   && n < lim
                x0 = x
                x = x0 - exp(lp - logpdf(d,x0) + log1mexp(min(logccdf(d,x0)-lp,0)))
                n += 1
            end
        end
        return x
    elseif lp == -Inf
        return T(maximum(d))
    elseif lp == 0
        return T(minimum(d))
    else
        return T(NaN)
    end
end