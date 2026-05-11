
abstract type AbstractMvStable <: ContinuousMultivariateDistribution end


struct EllipticStable{T<:Real,Cov<:AbstractPDMat,Mean<:AbstractVector} <: AbstractMvStable
    α::T
    Σ::Cov
    μ::Mean
    function EllipticStable(α::T, Σ::AbstractPDMat{T}, μ::AbstractVector{T}) where {T<:Real}
        0 < α <= 2 || throw(ArgumentError("Parameter α must fulfill 0 < α ≤ 2."))
        size(Σ, 1) == size(μ, 1) || throw(DimensionMismatch("The dimensions of μ and Σ are inconsistent."))
        new{T,typeof(Σ), typeof(μ)}(α, Σ, μ)
    end
end

function EllipticStable(α::T, Σ::AbstractPDMat{W}, μ::AbstractVector{S}) where {T<:Real, S<:Real, W<:Real}
    R = promote_type(T, S, W)
    EllipticStable(R(α), convert(AbstractArray{R}, Σ), convert(AbstractArray{R}, μ))
end

EllipticStable(α, Σ::AbstractMatrix{<:Real}, μ::AbstractVector{<:Real}) = 
    EllipticStable(α, PDMat(Σ), μ)
EllipticStable(α, Σ::Diagonal{<:Real}, μ::AbstractVector{<:Real}) = 
    EllipticStable(α, PDiagMat(Σ.diag), μ)
EllipticStable(α, Σ::Union{Symmetric{<:Real,<:Diagonal{<:Real}},Hermitian{<:Real,<:Diagonal{<:Real}}}, μ::AbstractVector{<:Real} ) = 
    EllipticStable(α, PDiagMat(Σ.data.diag), μ)
EllipticStable(α, Σ::UniformScaling{<:Real}, μ::AbstractVector{<:Real}) =
    EllipticStable(α, ScalMat(length(μ), Σ.λ), μ)
EllipticStable(α, Σ::Diagonal{<:Real,<:FillArrays.AbstractFill{<:Real,1}}, μ::AbstractVector{<:Real}) = 
    EllipticStable(α, ScalMat(size(Σ, 1), FillArrays.getindex_value(Σ.diag)), μ)

EllipticStable(α, Σ::AbstractMatrix{<:Real}) = 
    EllipticStable(α, Σ, Zeros{eltype(Σ)}(size(Σ, 1))) # default μ = [0, 0, ...]

insupport(d::EllipticStable, x::AbstractVector) =
    length(d) == length(x) && all(isfinite, x)

minimum(d::EllipticStable) = fill(typemin(eltype(d)), length(d))
maximum(d::EllipticStable) = fill(typemax(eltype(d)), length(d))


#### Parameters

length(d::EllipticStable) = length(d.μ)
params(d::EllipticStable) = (d.α, d.Σ, d.μ)
location(d::EllipticStable) = d.μ
mode(d::EllipticStable) = d.μ
Base.eltype(::Type{<:EllipticStable{T}}) where {T} = T

#### Conversion
convert(::Type{EllipticStable{T}}, d::EllipticStable) where T<:Real =
    EllipticStable(convert(T, d.α), convert(AbstractArray{T}, d.Σ), convert(AbstractArray{T}, d.μ))

Base.convert(::Type{EllipticStable{T}}, d::EllipticStable{T}) where {T<:Real} = d

convert(::Type{EllipticStable{T}}, α::Real, Σ::AbstractPDMat, μ::AbstractVector) where T<:Real =
    EllipticStable(convert(T, α), convert(AbstractArray{T}, Σ), convert(AbstractArray{T}, μ))



#### Affine transformations

Base.:+(d::EllipticStable, c::AbstractVector) = EllipticStable(d.α, d.Σ, μ + c)
Base.:+(c::AbstractVector, d::EllipticStable) = d + c
Base.:-(d::EllipticStable, c::AbstractVector) = EllipticStable(d.α, d.Σ, d.μ - c)

Base.:*(B::AbstractMatrix, d::EllipticStable) = EllipticStable(d.α, X_A_Xt(d.Σ, B), B * d.μ)

LinearAlgebra.dot(b::AbstractVector, d::EllipticStable) = Stable(d.α, 0, √quad(d.Σ, b), d.μ ⋅ b) # ?

LinearAlgebra.dot(d::EllipticStable, b::AbstractVector) = b ⋅ d

#### Evaluation

cf(d::EllipticStable, t::AbstractVector{<:Real}) = exp(-(quad(d.Σ, t))^(d.α/2) + im*d.μ ⋅ t)
cf(d::EllipticStable, t::AbstractMatrix{<:Real}) = exp.(-quad(d.Σ, t) .^(d.α/2) .+ im*t' * d.μ)


function _rand!(rng::AbstractRNG, d::EllipticStable, x::VecOrMat)
    unwhiten!(d.Σ, randn!(rng, x))
    x .= √rand(rng, Stable(d.α/2, 1, 2(cospi(d.α/4))^(2/d.α), 0)) .* x .+ d.μ
end

function _rand!(rng::AbstractRNG, d::EllipticStable, x::AbstractVector) 
    for i in eachindex(x)   
        x[i] = randn(rng, eltype(x))
    end
    unwhiten!(d.Σ, x)
    x .= √rand(rng, Stable(d.α/2, 1, 2(cospi(d.α/4))^(2/d.α), 0)) .* x .+ d.μ
end