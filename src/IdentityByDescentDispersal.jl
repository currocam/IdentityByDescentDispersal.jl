module IdentityByDescentDispersal
using BesselK: adbesselk
using QuadGK: quadgk
using DataFrames: DataFrame, transform
using Distributions: Poisson, logpdf
import Tables
import DataAPI
export safe_adbesselk,
    expected_ibd_blocks_constant_density,
    expected_ibd_blocks_power_density,
    expected_ibd_blocks_custom,
    preprocess_dataset,
    default_ibd_bins,
    composite_loglikelihood_constant_density,
    composite_loglikelihood_power_density,
    composite_loglikelihood_custom,
    age_density_ibd_blocks_custom

"""
    x = safe_adbesselk(1, 1.0)
A wrapper around the `adbesselk` to compute the modified Bessel function of the second kind of degree that returns a NaN if computation fails.
"""
@inline function safe_adbesselk(v, x)
    try
        return adbesselk(v, x)
    catch
        return NaN
    end
end

"""
    probability_coalescence(t::Real, r::Real, De::Function, sigma::Real)
Computes the probability ``\\phi(t)`` that two homologous loci coalesce ``t`` generations ago according to

```math
\\phi(t) = \\frac{1}{2D_e(t)} \\frac{1}{4\\pi t \\sigma^2} \\exp(\\frac{-r^2}{4 t \\sigma^2})
```

Ringbauer, H., Coop, G. & Barton, N.H. Genetics 205, 1335–1351 (2017).
"""
function probability_coalescence(t::Real, r::Real, De::Function, sigma::Real)
    probability_coalescence(t, r, De(t), sigma)
end

function probability_coalescence(t::Real, r::Real, De::Real, sigma::Real)
    if (r < 0 || De < 0 || sigma <= 0 || t < 0)
        throw(ArgumentError("Invalid input: r=$r, De=$De, sigma=$sigma, t=$t"))
    end
    if iszero(t)
        return zero(eltype(t))
    end
    sigma² = sigma^2
    denom = 4 * π * t * sigma²
    exponent = -r^2 / (4 * t * sigma²)
    return (1 / (2 * De)) * (1 / denom) * exp(exponent)
end

# Helper function from Appendix B
function nL(r::Real, D::Real, beta::Real, sigma::Real, L::Real)
    term1 = 2^(-3beta / 2 - 3)
    term2 = 1 / (π * D * sigma^2)
    term3 = (r / (sqrt(L) * sigma))^(2 + beta)
    term4 = safe_adbesselk(2 + beta, sqrt(2L) * r / sigma)
    term1 * term2 * term3 * term4
end

"""
    expected_ibd_blocks_constant_density(r::Real, D::Real, sigma::Real, L::Real, G::Real, chromosomal_edges::Bool=true, diploid::Bool=true) -> Real

Computes the expected density of identity-by-descent (IBD) blocks of length `L` for a model with constant population density.
This function returns the expected number of IBD blocks per pair of individuals and per unit of block length.

```math
\\mathbb{E}[N_L] = \\int_0^\\infty \\mathbb{E}[N_L^t] \\,dt =
\\frac{G}{8\\pi D \\sigma^2}
\\left(\\frac{r}{\\sqrt{L} \\sigma}\\right)^2
K_2\\left(\\frac{\\sqrt{2L} \\, r}{\\sigma}\\right)
```

where:
- `r` is the geographic distance between samples,
- `D` is the effective population density (diploid individuals per unit area),
- `sigma` is the root mean square dispersal distance per generation,
- `L` is the length of the IBD block (in Morgans),
- `G` is the total map length of the genome (in Morgans),
- `K₂` is the modified Bessel function of the second kind of order 2.

If `chromosomal_edges` is `true` (the default), we account for chromosomal edge effects. If `diploid` is `true`, we multiply by a factor of 4 to account for the fact that each individual has two copies of each chromosome. For more details, see Appendix B.

There is a function overload that accepts a vector of `G` values and returns the aggregated expected number of IBD blocks.

Reference:
Ringbauer, H., Coop, G., & Barton, N. H. (2017). Genetics, 205(3), 1335–1351.
"""
function expected_ibd_blocks_constant_density(
    r::Real,
    D::Real,
    sigma::Real,
    L::Real,
    G::Real,
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    # Input validation
    if r < 0 || D ≤ 0 || sigma ≤ 0 || L ≤ 0 || G ≤ 0
        throw(ArgumentError("All parameters must be positive"))
    end
    if chromosomal_edges
        result = (G - L) * nL(r, D, 0, sigma, L) + nL(r, D, -1, sigma, L)
    else
        term1 = G / (8π * D * sigma^2)
        term2 = (r / (sqrt(L) * sigma))^2
        term3 = safe_adbesselk(2, sqrt(2L) * r / sigma)
        result = term1 * term2 * term3
    end
    diploid ? result * 4 : result
end


function expected_ibd_blocks_constant_density(
    r::Real,
    D::Real,
    sigma::Real,
    L::Real,
    contig_lengths::AbstractArray{<:Real},
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    sum(
        expected_ibd_blocks_constant_density(
            r,
            D,
            sigma,
            L,
            contig_lengths[i],
            chromosomal_edges,
            diploid,
        ) for i in eachindex(contig_lengths)
    )
end

"""
    expected_ibd_blocks_power_density(r::Real, D::Real, beta::Real, sigma::Real, L::Real, G::Real, chromosomal_edges::Bool = true, diploid::Bool = true)

Computes the expected density of identity-by-descent (IBD) blocks of length `L` for a model with a power population density in the form of ``D(t) = D_0t^{-\beta}``
This function returns the expected number of IBD blocks per pair of individuals and per unit of block length.
```math
\\mathbb{E}[N_L] = \\int_0^\\infty \\mathbb{E}[N_L^t] \\,dt =
2^{\\frac{-3\\beta}{2}-3}
\\frac{G}{\\pi D \\sigma^2}
\\left(\\frac{r}{\\sqrt{L} \\sigma}\\right)^{(2+\\beta)}
K_{2+\\beta}\\left(\\frac{\\sqrt{2L} \\, r}{\\sigma}\\right)
```

where:
- `r` is the geographic distance between samples,
- `D` is the effective population density (diploid individuals per unit area),
- `beta` is the power of the density function,
- `sigma` is the root mean square dispersal distance per generation,
- `L` is the length of the IBD block (in Morgans),
- `G` is the total map length of the genome (in Morgans),
- `K₂` is the modified Bessel function of the second kind of order 2.

If `chromosomal_edges` is `true` (the default), we account for chromosomal edge effects. If `diploid` is `true`, we multiply by a factor of 4 to account for the fact that each individual has two copies of each chromosome. For more details, see Appendix B.

Reference:
Ringbauer, H., Coop, G., & Barton, N. H. (2017). Genetics, 205(3), 1335–1351.
"""
function expected_ibd_blocks_power_density(
    r::Real,
    D::Real,
    beta::Real,
    sigma::Real,
    L::Real,
    G::Real,
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    # Input validation
    if r < 0 || D ≤ 0 || sigma ≤ 0 || L ≤ 0 || G ≤ 0
        throw(ArgumentError("All parameters except beta must be positive"))
    end
    if chromosomal_edges
        result = (G - L) * nL(r, D, beta, sigma, L) + nL(r, D, beta - 1, sigma, L)
    else
        term1 = 2^(-3beta / 2 - 3)
        term2 = G / (π * D * sigma^2)
        term3 = (r / (sqrt(L) * sigma))^(2 + beta)
        term4 = safe_adbesselk(2 + beta, sqrt(2L) * r / sigma)
        result = term1 * term2 * term3 * term4
    end
    diploid ? result * 4 : result
end

function expected_ibd_blocks_power_density(
    r::Real,
    D::Real,
    beta::Real,
    sigma::Real,
    L::Real,
    contig_lengths::AbstractArray{<:Real},
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    sum(
        expected_ibd_blocks_power_density(
            r,
            D,
            beta,
            sigma,
            L,
            contig_lengths[i],
            chromosomal_edges,
            diploid,
        ) for i in eachindex(contig_lengths)
    )
end

"""
    expected_ibd_blocks_custom(r::Real, De::Function, parameters::AbstractArray, sigma::Real, L::Real, G::Real, chromosomal_edges::Bool = true, diploid::Bool = true)
Computes the expected density of identity-by-descent (IBD) blocks of length `L` for a model where the effective population density is given by a custom function `De(t, parameters)`.
This function returns the expected number of IBD blocks per pair of individuals and per unit of block length.

```math
\\mathbb{E}[N_L] = \\int_0^\\infty \\mathbb{E}[N_L^t] \\,dt = \\int_0^\\infty G  4t^2 \\exp(-2Lt) \\cdot \\Phi(t) \\,dt
```

The integral is computed numerically using Gaussian-Legendre quadrature rules with the `QuadGK` package.

where:
- `r` is the geographic distance between samples,
- `De` is  a user-defined function that takes time `t` and a `parameters` and returns the effective population density at time `t`.
- `parameters` is a user-defined array of parameters that the function `De` depends on.
- `sigma` is the root mean square dispersal distance per generation,
- `L` is the length of the IBD block (in Morgans),
- `G` is the total map length of the genome (in Morgans),
- `K₂` is the modified Bessel function of the second kind of order 2.

If `chromosomal_edges` is `true` (the default), we account for chromosomal edge effects. If `diploid` is `true`, we multiply by a factor of 4 to account for the fact that each individual has two copies of each chromosome. For more details, see Appendix B.

Reference:
Ringbauer, H., Coop, G., & Barton, N. H. (2017). Genetics, 205(3), 1335–1351.
"""
function expected_ibd_blocks_custom(
    r::Real,
    De::Function,
    parameters::AbstractArray,
    sigma::Real,
    L::Real,
    G::Real,
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    # Input validation
    if r < 0 || sigma ≤ 0 || L ≤ 0 || G ≤ 0
        throw(ArgumentError("All provided parameters must be positive"))
    end
    fn1(t) =
        probability_coalescence(t, r, De(t, parameters), sigma) *
        (4t * exp(-2L * t) + (G - L) * 4 * t^2 * exp(-2L * t))
    fn2(t) =
        probability_coalescence(t, r, De(t, parameters), sigma) * G * 4 * t^2 * exp(-2L * t)
    result = NaN
    try
        if chromosomal_edges
            result = quadgk(fn1, 0, Inf)[1]
        else
            result = quadgk(fn2, 0, Inf)[1]
        end
    catch e
        @warn "Integration failed with error: $e"
    end
    diploid ? result * 4 : result
end

function expected_ibd_blocks_custom(
    r::Real,
    De::Function,
    parameters::AbstractArray,
    sigma::Real,
    L::Real,
    contig_lengths::AbstractArray{<:Real},
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    sum(
        expected_ibd_blocks_custom(
            r,
            De,
            parameters,
            sigma,
            L,
            contig_lengths[i],
            chromosomal_edges,
            diploid,
        ) for i in eachindex(contig_lengths)
    )
end

"""
    default_ibd_bins()
Returns the default bins used in Ringbauer et. al. for identity-by-descent (IBD) block analysis, as well as the minimum length of IBD blocks to consider in Morgans.

Ringbauer, H., Coop, G. & Barton, N.H. Genetics 205, 1335–1351 (2017).
"""
function default_ibd_bins()
    minimum_threshold = 0.04 # 4cM
    maximum_threshold = 0.20 # 20cM
    step = 0.001 # 0.1 cM
    collect((minimum_threshold+step):step:maximum_threshold), minimum_threshold
end

"""
    preprocess_dataset(ibd_blocks::DataFrame, dist_matrix::AbstractMatrix, bins::AbstractVector, min_length::Real)

Preprocesses the input data for identity-by-descent (IBD) block analysis.

- `ibd_blocks`: DataFrame containing IBD blocks with columns `ID1`, `ID2`, and `span`. The `ID1` and `ID2` columns should contain the IDs of the individuals involved in the IBD block, and the `span` column should contain the length of the IBD block in Morgans.
- `individual_distances`: DataFrame containing distances between individuals with columns `ID1`, `ID2`, and `distance`. The `ID1` and `ID2` columns should contain the IDs of the individuals involved in the distance. Notice that the units of the estimated density and dispersal rate will match the units of the distances provided.
- `bins`: A vector of right bins for the IBD blocks.
- `min_length`: Minimum length of IBD blocks to consider in Morgans.

It returns a DataFrame in "long" format with the following columns:
  - `DISTANCE`: The pairwise distance between individuals.
  - `IBD_LEFT`: The left bound of the IBD length bin.
  - `IBD_RIGHT`: The right bound of the IBD length bin.
  - `IBD_MID`: The center of the IBD length bin.
  - `NR_PAIRS`: The number of unique individual pairs within that distance.
  - `COUNT`: The number of IBD blocks observed in the corresponding bin.
  - `DISTANCE_INDEX`: The index of the distance bin.
  - `IBD_INDEX`: The index of the IBD length bin.
"""
function preprocess_dataset(
    ibd_blocks::DataFrame,
    individual_distances::DataFrame,
    bins::AbstractVector,
    min_length::Real,
)
    # Normalize unordered pairs for matching: sort (ID1, ID2)
    function normalize_pair(id1, id2)
        id1 < id2 ? (id1, id2) : (id2, id1)
    end

    ibd_blocks =
        transform(ibd_blocks, [:ID1, :ID2] => Tables.ByRow(normalize_pair) => [:ID1, :ID2])
    individual_distances = transform(
        individual_distances,
        [:ID1, :ID2] => Tables.ByRow(normalize_pair) => [:ID1, :ID2],
    )

    # Result DataFrame
    result = DataFrame(
        DISTANCE = Float64[],
        IBD_LEFT = Float64[],
        IBD_RIGHT = Float64[],
        IBD_MID = Float64[],
        NR_PAIRS = Int[],
        COUNT = Int[],
        DISTANCE_INDEX = Int[],
        BIN_INDEX = Int[],
    )

    # Iterate over distances
    left_bins = [min_length; bins]
    bin_widths = bins .- bins
    for (dist_i, distance) in enumerate(unique(individual_distances.distance))
        # Get all pairs at this distance
        pairs_at_distance =
            individual_distances[individual_distances.distance .== distance, [:ID1, :ID2]]
        pair_set = Set(Tuple.(eachrow(pairs_at_distance)))

        # Filter IBD blocks for those pairs
        ibd_filtered = filter(row -> (row.ID1, row.ID2) ∈ pair_set, ibd_blocks)

        # Count unique pairs for this distance
        nr_pairs = length(pair_set)

        # Iterate over bins
        left_bin = min_length
        for (ibd_i, right_bin) in enumerate(bins)
            midpoint = (left_bin + right_bin) / 2
            # Select blocks in this bin
            in_bin = (ibd_filtered.span .>= left_bin) .& (ibd_filtered.span .< right_bin)
            ibd_in_bin = ibd_filtered[in_bin, :]
            count = DataAPI.nrow(ibd_in_bin)

            push!(
                result,
                (distance, left_bin, right_bin, midpoint, nr_pairs, count, dist_i, ibd_i),
            )
            left_bin = right_bin
        end
    end
    result
end


"""
    composite_loglikelihood_constant_density(D::Real, sigma::Real, df::DataFrame, contig_lengths::AbstractArray{<:Real}, chromosomal_edges::Bool=true, diploid::Bool=true) -> Real

Computes the composite log-likelihood of the observed IBD blocks under a model with constant population density.
- `D`: Effective population density (diploid individuals per unit area).
- `sigma`: Root mean square dispersal distance per generation.
- `df`: DataFrame containing the observed IBD blocks in the format returned by `preprocess_dataset`.
- `contig_lengths`: Array of contig lengths in Morgans.

Optionally:
- `chromosomal_edges`: Whether to account for chromosomal edge effects.
- `diploid`: Whether to account for diploidy.

"""
function composite_loglikelihood_constant_density(
    D::Real,
    sigma::Real,
    df::DataFrame,
    contig_lengths::AbstractArray{<:Real},
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    # Input validation
    if D ≤ 0 || sigma ≤ 0 || isempty(contig_lengths)
        return -Inf
    end

    # Iterate across each row computing the expected rate and updating the log-likelihood
    log_likelihood = 0.0
    for (i, row) in enumerate(eachrow(df))
        r, L, nr_pairs, count = row.DISTANCE, row.IBD_LEFT, row.NR_PAIRS, row.COUNT
        ΔL = row.IBD_RIGHT - row.IBD_LEFT
        λ =
            expected_ibd_blocks_constant_density(
                r,
                D,
                sigma,
                L,
                contig_lengths,
                chromosomal_edges,
                diploid,
            ) *
            ΔL *
            nr_pairs
        if λ < 0 || isnan(λ)
            return -Inf
        end
        log_likelihood += logpdf(Poisson(λ), count)
    end
    log_likelihood
end

"""
    composite_loglikelihood_power_density(D::Real, beta::Real, sigma::Real, df::DataFrame, contig_lengths::AbstractArray{<:Real}, chromosomal_edges::Bool=true, diploid::Bool=true) -> Real

Computes the composite log-likelihood of the observed IBD blocks under a model with constant population density.
- `D`: Effective population density (diploid individuals per unit area).
- `beta` is the power of the density function,
- `sigma`: Root mean square dispersal distance per generation.
- `df`: DataFrame containing the observed IBD blocks in the format returned by `preprocess_dataset`.
- `contig_lengths`: Array of contig lengths in Morgans.

Optionally:
- `chromosomal_edges`: Whether to account for chromosomal edge effects.
- `diploid`: Whether to account for diploidy.

"""
function composite_loglikelihood_power_density(
    D::Real,
    beta::Real,
    sigma::Real,
    df::DataFrame,
    contig_lengths::AbstractArray{<:Real},
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    # Input validation
    if D ≤ 0 || sigma ≤ 0 || isempty(contig_lengths)
        return -Inf
    end

    # Iterate across each row computing the expected rate and updating the log-likelihood
    log_likelihood = 0.0
    for (i, row) in enumerate(eachrow(df))
        r, L, nr_pairs, count = row.DISTANCE, row.IBD_LEFT, row.NR_PAIRS, row.COUNT
        ΔL = row.IBD_RIGHT - row.IBD_LEFT
        λ =
            expected_ibd_blocks_power_density(
                r,
                D,
                beta,
                sigma,
                L,
                contig_lengths,
                chromosomal_edges,
                diploid,
            ) *
            ΔL *
            nr_pairs
        if λ < 0 || isnan(λ)
            return -Inf
        end
        log_likelihood += logpdf(Poisson(λ), count)
    end
    log_likelihood
end

"""
    composite_loglikelihood_custom(De::Function, parameters::AbstractArray, sigma::Real, df::DataFrame, contig_lengths::AbstractArray{<:Real}, chromosomal_edges::Bool=true, diploid::Bool=true) -> Real

Computes the composite log-likelihood of the observed IBD blocks under a model with constant population density.
- `De` is  a user-defined function that takes time `t` and `parameters` and returns the effective population density at time `t`.
- `parameters` is a user-defined array of parameters that the function `De` depends on.
- `sigma`: Root mean square dispersal distance per generation.
- `df`: DataFrame containing the observed IBD blocks in the format returned by `preprocess_dataset`.
- `contig_lengths`: Array of contig lengths in Morgans.

Optionally:
- `chromosomal_edges`: Whether to account for chromosomal edge effects.
- `diploid`: Whether to account for diploidy.

"""
function composite_loglikelihood_custom(
    De::Function,
    parameters::AbstractArray,
    sigma::Real,
    df::DataFrame,
    contig_lengths::AbstractArray{<:Real},
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    # Input validation
    if sigma ≤ 0 || isempty(contig_lengths)
        return -Inf
    end

    # Iterate across each row computing the expected rate and updating the log-likelihood
    log_likelihood = 0.0
    for (i, row) in enumerate(eachrow(df))
        r, L, nr_pairs, count = row.DISTANCE, row.IBD_LEFT, row.NR_PAIRS, row.COUNT
        ΔL = row.IBD_RIGHT - row.IBD_LEFT
        λ =
            expected_ibd_blocks_custom(
                r,
                De,
                parameters,
                sigma,
                L,
                contig_lengths,
                chromosomal_edges,
                diploid,
            ) *
            ΔL *
            nr_pairs
        if λ < 0 || isnan(λ)
            return -Inf
        end
        log_likelihood += logpdf(Poisson(λ), count)
    end
    log_likelihood
end

# For posterior predictive simulations

"""
    age_density_ibd_blocks_custom(t::Real, r::Real, De::Function, parameters::AbstractArray, sigma::Real, L::Real, G::Real, chromosomal_edges::Bool = true, diploid::Bool = true)
Computes the expected density of identity-by-descent (IBD) blocks of length `L` and age `t` for a model where the effective population density is given by a custom function `De(t, parameters)`.
This function returns the expected number of IBD blocks of age `t` per pair of individuals and per unit of block length.

```math
\\mathbb{E}[N_L^t] =  G  4t^2 \\exp(-2Lt) \\cdot \\Phi(t)
```


where:
- `t` is the age of the IBD block,
- `r` is the geographic distance between samples,
- `De` is a user-defined function that takes time `t` and a `parameters` and returns the effective population density at time `t`.
- `parameters` is a user-defined array of parameters that the function `De` depends on.
- `sigma` is the root mean square dispersal distance per generation,
- `L` is the length of the IBD block (in Morgans),
- `G` is the total map length of the genome (in Morgans),

If `chromosomal_edges` is `true` (the default), we account for chromosomal edge effects. If `diploid` is `true`, we multiply by a factor of 4 to account for the fact that each individual has two copies of each chromosome. For more details, see Appendix B.

Reference:
Ringbauer, H., Coop, G., & Barton, N. H. (2017). Genetics, 205(3), 1335–1351.
"""
function age_density_ibd_blocks_custom(
    t::Real,
    r::Real,
    De::Function,
    parameters::AbstractArray,
    sigma::Real,
    L::Real,
    G::Real,
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    # Input validation
    if r < 0 || sigma ≤ 0 || L ≤ 0 || G ≤ 0
        throw(ArgumentError("All provided parameters must be positive"))
    end
    fn(t) = begin
        phi = probability_coalescence(t, r, De(t, parameters), sigma)
        if chromosomal_edges
            phi * (4t * exp(-2L * t) + (G - L) * 4 * t^2 * exp(-2L * t))
        else
            phi * G * 4 * t^2 * exp(-2L * t)
        end
    end
    result = fn(t)
    diploid ? result * 4 : result
end

function age_density_ibd_blocks_custom(
    t::Real,
    r::Real,
    De::Function,
    parameters::AbstractArray,
    sigma::Real,
    L::Real,
    G::AbstractArray,
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    sum([
        age_density_ibd_blocks_custom(
            t,
            r,
            De,
            parameters,
            sigma,
            L,
            g,
            chromosomal_edges,
            diploid,
        ) for g in G
    ])
end

function age_density_ibd_blocks_custom(
    t::AbstractArray,
    r::Real,
    De::Function,
    parameters::AbstractArray,
    sigma::Real,
    L::Real,
    G::Union{Real,AbstractArray},
    chromosomal_edges::Bool = true,
    diploid::Bool = true,
)
    [
        age_density_ibd_blocks_custom(
            ti,
            r,
            De,
            parameters,
            sigma,
            L,
            G,
            chromosomal_edges,
            diploid,
        ) for ti in t
    ]
end


end
