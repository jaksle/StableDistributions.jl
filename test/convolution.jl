using Distributions
using Test



@testset "continuous univariate" begin
    @testset "Stable" begin
        @test_throws ArgumentError convolve(Stable(1),Stable(1.5))
        d1 = Stable(1.5, 1/2,2, 1)
        d2 = Stable(1.5, -1/2, 1, 1)
        d3 = @inferred(convolve(d1, d2))
        @test d3 isa Stable
        @test d3.α == 1.5
        @test d3.β ≈ (d1.β*d1.σ^d3.α + d2.β*d2.σ^d3.α) / (d1.σ^d3.α + d2.σ^d3.α)
        @test d3.σ ≈ (d1.σ^d3.α + d2.σ^d3.α)^(1/d3.α)
        @test d3.μ == 1 + 1
    end
end
