using IdentityByDescentDispersal
using Random, Test, Distributions, QuadGK


@testset "IdentityByDescentDispersal.jl" begin
    @testset "Testing pointwise coalescent probabilities" begin
        # A hard-coded example
        expected = (1 / (2 * 1)) * (1 / (4 * π * 1 * 1^2)) * exp(-1^2 / (4 * 1 * 1^2))
        computed = IdentityByDescentDispersal.ϕ(1, 1, t -> 1, 1)
        @test isapprox(computed, expected; atol = 1e-5)
        # Invariance tests:
        for _ = 1:50
            r = rand(Uniform(0, 1e5))
            σ = rand(Uniform(0, 1e5))
            de = t -> rand(Uniform(0, 1e5))
            @test iszero(IdentityByDescentDispersal.ϕ(0, r, de, σ))
            computed = IdentityByDescentDispersal.ϕ.(rand(Uniform(0, 1e5), 100), r, de, σ)
            @test all(computed .>= 0)
            @test all(computed .<= 1)
        end
    end
    @testset "Testing analytical expected IBD blocks versus numerical integration" begin
        for _ = 1:50
            r = rand(Uniform(0, 10))
            D = rand(Uniform(0, 10))
            σ = rand(Uniform(0, 10))
            L = rand(Uniform(0, 0.5))
            G = rand(Uniform(L, 1))
            # Constant population size (without chromosomal edge effects)
            analytical_value =
                IdentityByDescentDispersal.expected_ibd_blocks_constant_density(
                    r,
                    D,
                    σ,
                    L,
                    G,
                    false,
                    false,
                )
            fn = x -> D
            numerical_value = quadgk(
                t ->
                    IdentityByDescentDispersal.ϕ(t, r, fn, σ) * G * 4 * t^2 * exp(-2L * t),
                0,
                Inf,
            )[1]
            @test isapprox(numerical_value, analytical_value; atol = 1e-5)
            # Constant population size (with chromosomal edge effects)
            analytical_value =
                IdentityByDescentDispersal.expected_ibd_blocks_constant_density(
                    r,
                    D,
                    σ,
                    L,
                    G,
                    true,
                    false,
                )
            fn = x -> D
            numerical_value = quadgk(
                t ->
                    IdentityByDescentDispersal.ϕ(t, r, fn, σ) *
                    (4t * exp(-2L * t) + (G - L) * 4 * t^2 * exp(-2L * t)),
                0,
                Inf,
            )[1]
            @test isapprox(numerical_value, analytical_value; atol = 1e-5)

        end
        # Power density
        for _ = 1:50
            r = rand(Uniform(0, 10))
            D = rand(Uniform(0, 10))
            σ = rand(Uniform(0, 10))
            L = rand(Uniform(0, 0.5))
            G = rand(Uniform(L, 1))
            β = rand(Normal(0, 1))
            #  without chromosomal edge effects
            analytical_value = IdentityByDescentDispersal.expected_ibd_blocks_power_density(
                r,
                D,
                β,
                σ,
                L,
                G,
                false,
                false,
            )
            fn = x -> D * x^(-β)
            numerical_value = quadgk(
                t ->
                    IdentityByDescentDispersal.ϕ(t, r, fn, σ) * G * 4 * t^2 * exp(-2L * t),
                0,
                Inf,
            )[1]
            @test isapprox(numerical_value, analytical_value; atol = 1e-5)
            # with chromosomal edge effects
            analytical_value = IdentityByDescentDispersal.expected_ibd_blocks_power_density(
                r,
                D,
                β,
                σ,
                L,
                G,
                true,
                false,
            )
            fn = x -> D * x^(-β)
            numerical_value = quadgk(
                t ->
                    IdentityByDescentDispersal.ϕ(t, r, fn, σ) *
                    (4t * exp(-2L * t) + (G - L) * 4 * t^2 * exp(-2L * t)),
                0,
                Inf,
            )[1]
            @test isapprox(numerical_value, analytical_value; atol = 1e-5)
        end
        # Constant density via user-defined function
        for _ = 1:50
            r = rand(Uniform(0, 10))
            D1 = rand(Uniform(0, 5))
            D2 = rand(Uniform(0, 5))
            σ = rand(Uniform(0, 10))
            L = rand(Uniform(0, 0.5))
            G = rand(Uniform(L, 1))
            fn(t, θ) = θ[1] + θ[2]
            # Without chromosomal edge effects
            computed = IdentityByDescentDispersal.expected_ibd_blocks_custom(
                r,
                fn,
                [D1, D2],
                σ,
                L,
                G,
                false,
                false,
            )
            expected = IdentityByDescentDispersal.expected_ibd_blocks_constant_density(
                r,
                D1 + D2,
                σ,
                L,
                G,
                false,
                false,
            )
            @test isapprox(computed, expected; atol = 1e-5)
            # With chromosomal edge effects
            computed = IdentityByDescentDispersal.expected_ibd_blocks_custom(
                r,
                fn,
                [D1, D2],
                σ,
                L,
                G,
                true,
                false,
            )
            expected = IdentityByDescentDispersal.expected_ibd_blocks_constant_density(
                r,
                D1 + D2,
                σ,
                L,
                G,
                true,
                false,
            )
            @test isapprox(computed, expected; atol = 1e-5)

        end
    end
end
