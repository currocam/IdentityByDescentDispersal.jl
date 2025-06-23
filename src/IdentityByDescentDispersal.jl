module IdentityByDescentDispersal
using BesselK: adbesselk
using QuadGK: quadgk
using DataFrames: DataFrame, transform
import Tables
import DataAPI
export safe_adbesselk,
    expected_ibd_blocks_constant_density,
    expected_ibd_blocks_power_density,
    expected_ibd_blocks_custom,
    preprocess_dataset

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
    if iszero(t)
        return zero(eltype(t))
    end
    sigma² = sigma^2
    denom = 4 * π * t * sigma²
    exponent = -r^2 / (4 * t * sigma²)
    return (1 / (2 * De(t))) * (1 / denom) * exp(exponent)
end

function probability_coalescence(t::Real, r::Real, De::Real, sigma::Real)
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

Computes the expected number of identity-by-descent (IBD) blocks of length `L` for a model with constant population density.

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
- `L` is the minimum length of the IBD block (in Morgans),
- `G` is the total map length of the genome (in Morgans),
- `K₂` is the modified Bessel function of the second kind of order 2.

If `chromosomal_edges` is `true` (the default), we account for chromosomal edge effects. If `diploid` is `true`, we multiply by a factor of 4 to account for the fact that each individual has two copies of each chromosome. For more details, see Appendix B.

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

"""
    expected_ibd_blocks_power_density(r::Real, D::Real, beta::Real, sigma::Real, L::Real, G::Real, chromosomal_edges::Bool = true, diploid::Bool = true)

Computes the expected number of identity-by-descent (IBD) blocks of length `L` for a model with a power population density in the form of ``D(t) = D_0t^{-beta}``

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
- `beta` is the growth rate of the population
- `sigma` is the root mean square dispersal distance per generation,
- `L` is the minimum length of the IBD block (in Morgans),
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


"""
    expected_ibd_blocks_custom(r::Real, De::Function, parameters::AbstractArray, sigma::Real, L::Real, G::Real, chromosomal_edges::Bool = true, diploid::Bool = true)
Computes the expected number of identity-by-descent (IBD) blocks of length `L` for a model where the effective population density is given by a custom function `De(t, parameters)`.

```math
\\mathbb{E}[N_L] = \\int_0^\\infty \\mathbb{E}[N_L^t] \\,dt = \\int_0^\\infty \\mathbb{E}[N_L^t] \\,dt  G  4t^2 \\exp(-2Lt) \\cdot \\Phi(t) \\,dt
```

The integral is computed numerically using Gaussian-Legendre quadrature rules with the `QuadGK` package.

where:
- `r` is the geographic distance between samples,
- `De` is  a user defined function that takes time `t` and a `parameters` and returns the effective population density at time `t`.
- `parameters` is a user defined array of parameters that the function `De` depends on.
- `sigma` is the root mean square dispersal distance per generation,
- `L` is the minimum length of the IBD block (in Morgans),
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
    if chromosomal_edges
        result = quadgk(fn1, 0, Inf)[1]
    else
        result = quadgk(fn2, 0, Inf)[1]
    end
    diploid ? result * 4 : result
end

"""
    preprocess_dataset(ibd_blocks::DataFrame, dist_matrix::AbstractMatrix, bins::AbstractVector, min_length::Real)

Preprocesses the input data for identity-by-descent (IBD) block analysis.

- `ibd_blocks`: DataFrame containing IBD blocks with columns `ID1`, `ID2`, and `span`. The `ID1` and `ID2` columns should contain the IDs of the individuals involved in the IBD block, and the `span` column should contain the length of the IBD block in Morgans.
- `individual_distances`: DataFrame containing distances between individuals with columns `ID1`, `ID2`, and `distance`. The `ID1` and `ID2` columns should contain the IDs of the individuals involved in the distance. Notice that the units of the estimated density and dispersal rate will match the units of the distances provided.
- `bins`: A vector of right bins for the IBD blocks.
- `min_length`: Minimum length of IBD blocks to consider in Morgans.

It returns a DataFrame in "long" format with the following columns:
  - `DISTANCE`: The pairwise distance bin between individuals.
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
    IBDDispersalDataset
A container for storing preprocessed data required for the inference of dispersal rates


- `contig_lengths::Vector{Float64}`: The lengths of each contig (in Morgans) used to normalize IBD span counts, if applicable.
- `bin_widths::Vector{Float64}`: The widths of the IBD length bins, computed from `bin_edges`.

"""
struct IBDDispersalDataset
    ibd_summary::DataFrame
    contig_lengths::Vector{Float64}
    bin_widths::Vector{Float64}
end

end
