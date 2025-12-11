
convert(::Type{Stable}, d::Normal) = Stable(2, 0, d.σ/√2, d.μ)
convert(::Type{Stable}, d::Cauchy) = Stable(1, 0, d.σ, d.μ)
convert(::Type{Stable}, d::Levy) = Stable(1/2, 1, d.σ, d.μ)

