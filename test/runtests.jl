using IdentityByDescentDispersal
using Random, Test, Distributions, QuadGK, DataFrames


@testset "IdentityByDescentDispersal.jl" begin
    @testset "Testing pointwise coalescent probabilities" begin
        # A hard-coded example
        expected = (1 / (2 * 1)) * (1 / (4 * π * 1 * 1^2)) * exp(-1^2 / (4 * 1 * 1^2))
        computed = IdentityByDescentDispersal.probability_coalescence(1, 1, t -> 1, 1)
        @test isapprox(computed, expected; atol = 1e-5)
        # Invariance tests:
        for _ = 1:50
            r = rand(Uniform(0, 1e5))
            σ = rand(Uniform(0, 1e5))
            de = t -> rand(Uniform(0, 1e5))
            @test iszero(IdentityByDescentDispersal.probability_coalescence(0, r, de, σ))
            computed =
                IdentityByDescentDispersal.probability_coalescence.(
                    rand(Uniform(0, 1e5), 100),
                    r,
                    de,
                    σ,
                )
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
                    IdentityByDescentDispersal.probability_coalescence(t, r, fn, σ) *
                    G *
                    4 *
                    t^2 *
                    exp(-2L * t),
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
                    IdentityByDescentDispersal.probability_coalescence(t, r, fn, σ) *
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
                    IdentityByDescentDispersal.probability_coalescence(t, r, fn, σ) *
                    G *
                    4 *
                    t^2 *
                    exp(-2L * t),
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
                    IdentityByDescentDispersal.probability_coalescence(t, r, fn, σ) *
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
    @testset "Testing data preparation" begin
        # We simply compare against a hardcoded example
        # Input data
        ibd_blocks = DataFrame(
            ID1 = ["A", "B", "A", "C", "B"],
            ID2 = ["B", "A", "C", "A", "C"],
            span = [0.006, 0.012, 0.021, 0.010, 0.008],
        )
        individual_distances = DataFrame(
            ID1 = ["A", "A", "B"],
            ID2 = ["B", "C", "C"],
            distance = [10.0, 20.0, 10.0],
        )
        bins = [0.01, 0.02]
        min_length = 0.005
        # Run the function
        computed = preprocess_dataset(ibd_blocks, individual_distances, bins, min_length)
        # Harcoded DataFrame
        expected = DataFrame(
            DISTANCE = vcat(repeat([10.0], 2), repeat([20.0], 2)),
            IBD_LEFT = repeat([0.005, 0.01], 2),
            IBD_RIGHT = repeat([0.01, 0.02], 2),
            IBD_MID = repeat([0.0075, 0.015], 2),
            NR_PAIRS = Int[2, 2, 1, 1],
            COUNT = Int[2, 1, 0, 1],
            DISTANCE_INDEX = Int[1, 1, 2, 2],
            BIN_INDEX = Int[1, 2, 1, 2],
        )
        @test computed == expected
    end
    @testset "Test that composite likelihood function do not throw errors" begin
        n_obs = 10
        left_bins = sort(rand(Uniform(0, 0.10), n_obs))
        delta = rand(Uniform(0, 0.01), n_obs)
        right_bins = left_bins .+ delta
        n_contigs = 3
        contig_lengths = rand(Uniform(0.5, 1.5), n_contigs)
        df = DataFrame(
            DISTANCE = rand(Uniform(0, 10), n_obs),
            IBD_LEFT = left_bins,
            IBD_RIGHT = right_bins,
            NR_PAIRS = rand(1:50, n_obs),
            COUNT = rand(0:50, n_obs),
        )
        for _ = 1:50
            D = rand(Uniform(0, 5))
            σ = rand(Uniform(0, 10))
            beta = rand(Uniform(-1, 1))
            alpha = rand(Normal(0, 0.001))
            isnotnan(x) = !isnan(x)
            # Constant population scenario
            @test isnotnan(
                composite_loglikelihood_constant_density(D, σ, df, contig_lengths),
            )
            # Power density
            @test isnotnan(
                composite_loglikelihood_power_density(D, beta, σ, df, contig_lengths),
            )
            # Exponential density
            De(t, parameters) = max(parameters[1] * exp(-parameters[2] * t), 1e-5)
            @test isnotnan(
                composite_loglikelihood_custom(De, [D, alpha], σ, df, contig_lengths),
            )
        end
    end

    @testset "Testing posterior predictive" begin
        r = rand(Uniform(0, 20))
        De(t, params) = t * params[1]
        sigma = rand(Uniform(0.1, 10.0))
        L = rand(Uniform(1e-4, 1e-2))
        t1 = rand(Uniform(5.0, 15.0))
        t2 = rand(Uniform(5.0, 15.0))
        param1 = rand(Uniform(1.0, 20.0))
        G1 = rand(Uniform(0.5, 2.0))
        G2 = rand(Uniform(0.5, 2.0))
        # For t=t1 and G=G1
        val = age_density_ibd_blocks_custom(t1, r, De, [param1], sigma, L, G1)
        @test isa(val, Number)
        # For t=t1 and G=[G1, G2]
        val = age_density_ibd_blocks_custom(t1, r, De, [param1], sigma, L, [G1, G2])
        @test isa(val, Number)
        # For t=[t1, t2] and G=G1
        val = age_density_ibd_blocks_custom([t1, t2], r, De, [param1], sigma, L, G1)
        @test isa(val, Vector{Float64})
        # For t=[t1, t2] and G=[G1, G2]
        val = age_density_ibd_blocks_custom([t1, t2], r, De, [param1], sigma, L, [G1, G2])
        @test isa(val, Vector{Float64})
        # Operations should be vectorized through time
        res1 = [
            age_density_ibd_blocks_custom(t1, r, De, [param1], sigma, L, [G1, G2]),
            age_density_ibd_blocks_custom(t2, r, De, [param1], sigma, L, [G1, G2]),
        ]
        res2 = age_density_ibd_blocks_custom([t1, t2], r, De, [param1], sigma, L, [G1, G2])
        @test res1 ≈ res2
        # The density across multiple contigs should be summed
        res1 = sum([
            age_density_ibd_blocks_custom(t1, r, De, [param1], sigma, L, G1),
            age_density_ibd_blocks_custom(t1, r, De, [param1], sigma, L, G1),
        ])
        res2 = age_density_ibd_blocks_custom(t1, r, De, [param1], sigma, L, [G1, G1])
        @test sum(res1) ≈ sum(res2)
    end
end
