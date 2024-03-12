

function convolve(d1::Stable, d2::Stable)
    d1.α ≈ d2.α || throw(ArgumentError("$(d1.α) !≈ $(d2.α): α parameters must be approximately equal"))
    α, β₁, σ₁, μ₁ = params(d1)
    _, β₂, σ₂, μ₂ = params(d2)

    return Stable(
        α,                                   # α
        (β₁*σ₁^α + β₂*σ₂^α) / (σ₁^α + σ₂^α), # β
        (σ₁^α + σ₂^α)^(1/α),                 # σ
        μ₁ + μ₂                              # μ
    )
end
