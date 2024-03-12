using Test
using Distributions

@test convert(Stable, Normal(42., 2.)) ≈ Stable(2., 0., √2, 42.)
@test convert(Stable, Cauchy(7., 13.)) == Stable(1., 0., 13., 7.)
@test convert(Stable, Levy(3., 5.)) == Stable(0.5, 1., 5., 3.)
