# This section demonstrates the basic usage of the package. We refer to the `Inference and model evaluation` section
# of the documentation for an explanation of how to use this package to estimate parameters of demographic models from observed data.

using IdentityByDescentDispersal

# ## Expected density of identity-by-descent blocks

# Genetic distances between populations are often a function of geographic distance. This phenomenon is often referred to as isolation-by-distance.
# Isolation-by-distance is reflected in the patterns of identity-by-descent blocks (also known as IBD blocks).
# It is possible to derive the expected patterns of shared identity-by-descent blocks under a diffusion approximation and certain demographic models.
# We refer to the publication of Ringbauer et al. (2017) and the `Model overview` section of the documentation for more details.
#
# The main feature of this package is to compute the expected density of identity-by-descent blocks under different demographic scenarios.
# Next, we show how to compute the expected density of identity-by-descent blocks under a constant population density using the `expected_ibd_blocks_constant_density` function.

using Plots
using QuadGK

L = 0.01         # Block length threshold (in Morgans)
G = 1.0          # Genome length (in Morgans)
r_values = range(0.01, 25.0, length = 200);  # Distances

plot(
    r_values,
    expected_ibd_blocks_constant_density.(r_values, 1.0, 0.5, L, G),
    xlabel = "Geographic distance between individuals (r)",
    ylabel = "E[#IBD blocks per unit of block length and per pair]",
    title = "Constant effective density scenario",
    label = "D=1.0, σ=0.5",
)
plot!(
    r_values,
    expected_ibd_blocks_constant_density.(r_values, 2.0, 0.5, L, G),
    label = "D=2.0, σ=0.5",
)
plot!(
    r_values,
    expected_ibd_blocks_constant_density.(r_values, 1.0, 0.8, L, G),
    label = "D=1.0, σ=0.8",
)

# The figure above shows how the density of IBD blocks of length 1 centimorgan in a genome of 1 Morgan decays as a function of geographic distance
# for different demographic scenarios that vary in effective density and dispersal rate.
# The effective density is in units of number of diploid individuals per area unit (e.g., km²). The dispersal rate is defined as the square root
# of the average squared axial parent-offspring distance. The rate of the decay of identity-by-descent blocks is directly related
# to the time since the most recent common ancestor.
#
# Recall that the length of IBD blocks is a continuous random variable. Therefore, the expected number of IBD blocks of exactly length `L` is zero.
# What the `expected_ibd_blocks_constant_density` function returns is the expected number of IBD blocks of length `L` per unit of block
# length (in Morgans) and per pair. In order to compare such expectations with observed data, we often want to consider instead the
# expected number of blocks whose length falls in a given interval `[a, b]`. For a small enough interval $[L, L + \Delta L]$ we can approximate
# the expected number of blocks whose length falls in that interval as:
#
# $E[N_{[L, L + \Delta L]}] = \int_{L}^{L + \Delta L} E[N_L] dL \approx E[N_L] \Delta L$
#
# For example, the expected number of identity-by-descent blocks whose length falls in the interval `[1cM, 1.1cM]` shared by two diploid
# individuals that are 1 kM apart in a population with constant effective population density of `10/kM²`and dispersal rate of `0.5/kM` would be

expected_ibd_blocks_constant_density(1.0, 10.0, 0.5, L, G) * 0.001

# However, we can also compute the expected number of blocks for any bin size using numerical integration (and we do this when computing log-likelihoods internally when necessary):

quadgk(x -> expected_ibd_blocks_constant_density(x, 10.0, 0.5, L, G) * 0.001, 1.0, 1.1)[1]

# Ringbauer et al. (2017) derived analytical solutions for the expected density of identity-by-descent blocks for
# the family of functions of effective population density of the form $D_e(t) = D_1 t^{-\beta}$. A constant effective population density is a special
# case of this family with $\beta = 0$. For other values of $\beta$, the expected density of identity-by-descent blocks can be computed using the
# `expected_ibd_blocks_power_density` function.

# In addition to power density functions, we also support arbitrary effective density functions provided by the user via numerical integration.
# Next, we illustrate how to define a custom function by considering the popular choice of an exponential growth function.

D = 1.0          # Effective population density
σ = 1.0          # Dispersal rate
function De(t, θ)
    D₀, α = θ
    max(D₀ * exp(-α * t), eps())
end
plot(
    r_values,
    [expected_ibd_blocks_custom(r, De, [D, 0.0], σ, L, G) for r in r_values],
    xlabel = "Geographic distance between individuals (r)",
    ylabel = "E[#IBD blocks per unit of block length and per pair]",
    title = "Constant effective density scenario",
    label = "growth rate α = 0.0",
)
plot!(
    r_values,
    [expected_ibd_blocks_custom(r, De, [D, 0.005], σ, L, G) for r in r_values],
    label = "growth rate α = 0.005",
)
plot!(
    r_values,
    [expected_ibd_blocks_custom(r, De, [D, -0.005], σ, L, G) for r in r_values],
    label = "growth rate α = -0.005",
)

# It is possible to define more complex effective density functions. However, users should consider carefully whether such parameters
# are identifiable from identity-by-descent blocks if they aim to estimate them. As an example, we compare the one scenario with an
# oscillating effective density function with a constant effective density demography.

function De(t, θ)
    D₀, a, ω = θ
    D₀ * (1 + a * sin(ω * t))
end
D = 1.0          # Effective population density
σ = 1.0          # Dispersal rate
θ = [D, 0.5, 2π]  # Parameters for De(t): D₀, a, ω
t_values = range(0.0, 10.0, length = 200)   # Time for plotting D_e(t)
p1 = plot(
    t_values,
    [De(t, θ) for t in t_values],
    xlabel = "Time (generations ago)",
    ylabel = "Effective population density",
    label = "Oscillating effective density",
)
hline!(p1, [D], label = "Constant effective density")
p2 = plot(
    r_values,
    [expected_ibd_blocks_custom(r, De, θ, σ, L, G) for r in r_values],
    xlabel = "Geographic distance between individuals (r)",
    ylabel = "E[#IBD blocks per unit of block length and per pair]",
    label = "Oscillating effective density",
)
plot!(
    p2,
    r_values,
    [expected_ibd_blocks_constant_density(r, D, σ, L, G) for r in r_values],
    label = "Constant effective density",
)

plot(p1, p2, layout = (2, 1), size = (600, 800))


# ## Inference

# Ringbauer et. al (2017) proposed to do inference by assuming that the observed number of IBD blocks that a pair $r$ units apart that fall within a small bin $[L_i, L_i + \Delta L]$ follows a Poisson distribution with mean $E[N_{L_i}(r, \theta)] \Delta L$.

# Therefore, the log probability of observing `Y` identity-by-descent blocks whose length fall in the bin $[L_i, L_i + \Delta L]$ from a pair of individuals that are $r$ units apart can be calculated simply as:

using Distributions
λ = expected_ibd_blocks_constant_density(
    0.2, # Distance,
    0.5, # Dispersal rate
    2.0, # Effective density
    0.01, # Block length threshold (in Morgans)
    1.0, # Genome length (in Morgans)
)
Y = 5 # Observed number of IBD blocks
logpdf(Poisson(λ), Y)

# The inference scheme is based on a composite-likelihood that treats each bin from different pairs of individuals as independent. For computational reasons, it is also faster to aggregate observations that correspond to the same distance and bin.

# The input data for this sort of analysis is:

# - A `DataFrame` containing the length of identity-by-descent blocks shared across different individuals.
# - A `DataFrame` containing distances across pairs of individuals.
# - A list of contig lengths to properly take into account chromosomal edges (in Morgans).

using DataFrames
ibd_blocks = DataFrame(
    ID1 = ["A", "B", "A", "C", "B"],
    ID2 = ["B", "A", "C", "A", "C"],
    span = [0.005, 0.012, 0.21, 0.10, 0.08], # Morgans
)
individual_distances = DataFrame(
    ID1 = ["A", "A", "B"],
    ID2 = ["B", "C", "C"],
    distance = [10.0, 20.0, 10.0], # (e.g., in kilometers)
)
contig_lengths = [1.0, 1.5]; # Contig lengths (in Morgans)

# Preprocessing the data requires combining the different sources of information and binning the identity-by-descent blocks.
# For this purpose, we provide a helper function `preprocess_dataset` that returns a `DataFrame` in a "long" format.
# Of course, you may specify your own bins, but here we use the same bins used by Ringbauer et al. (2017), which can be accessed via the `default_ibd_bins` function.
bins, min_length = default_ibd_bins()
df = preprocess_dataset(ibd_blocks, individual_distances, bins, min_length);

# In most scenarios, distances between diploid individuals are not available, but between sampling sites. The `preprocess_dataset` function handles this scenario appropriately
# by aggregating observations from pairs that have the exact same distance. However, if distances are available at the individual level, you may consider binning them to reduce
# the runtime. You might find an example of this in the [simulations](https://github.com/currocam/IdentityByDescentDispersal.jl/tree/main/simulations) subdirectory.

# We can now use the df `DataFrame`to compute composite log likelihoods of different parameters using the family of `composite_loglikelihood_*` functions.
# The second major feature of this package is that functions are all compatible with automatic differentiation. Analytical solutions depend on a modified Bessel function of the second kind, which is not widely implemented in a way that allows for automatic differentiation .
# Here, we rely on the `BesselK.jl` package for doing this. In addition, we rely on the system of `QuadGK.jl` for numerical integration compatible with automatic differentiation. This allows us to perform gradient-based optimization and interact with existing software.
#
# We recommend using this package together with the `Turing.jl`, as it provides a convenient interface for Bayesian and frequentist inference.
#
# ### Maximum likelihood estimation

# For example, we can estimate the parameters of a constant density model using maximum likelihood estimation:

using Turing, StatsBase, StatsPlots
@model function constant_density(df, contig_lengths)
    D ~ Uniform(0, 100)
    σ ~ Uniform(0, 20)
    Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
end

# Generate a MLE estimate.
mle_estimate = maximum_likelihood(constant_density(df, contig_lengths))
DataFrame(coeftable(mle_estimate)) # computed from the Fisher information matrix

# ### Bayesian inference
#
# Alternatively, we can do a standard Bayesian inference with any of the available inference algorithms such as NUTS. Here, we fit a power-density model using NUTS:

@model function power_density(df, contig_lengths)
    D ~ Truncated(Normal(100, 20), 0, Inf)
    σ ~ Truncated(Normal(1, 0.1), 0, Inf)
    β ~ Normal(0, 0.5)
    Turing.@addlogprob! composite_loglikelihood_power_density(D, β, σ, df, contig_lengths)
end
m = power_density(df, contig_lengths)
# Sample 4 chains in a serial fashion.
chains = sample(m, NUTS(), MCMCSerial(), 1000, 4)
plot(chains)

# More complex models can be built with cautious consideration of identifiability. For example, we can fit a piecewise exponential density model using a custom function and numerical integration.

function piecewise_D(t, parameters)
    De1, De2, alpha, t0 = parameters
    t <= t0 ? De1 * exp(-alpha * t) : De2
end
@model function exponential_density(df, contig_lengths)
    De1 ~ Truncated(Normal(1000, 100), 0, Inf)
    De2 ~ Truncated(Normal(1000, 100), 0, Inf)
    α ~ Normal(0, 0.05)
    t0 ~ Truncated(Normal(500, 100), 0, Inf)
    σ ~ Truncated(Normal(1, 0.1), 0, Inf)
    theta = [De1, De2, α, t0]
    Turing.@addlogprob! composite_loglikelihood_custom(
        piecewise_D,
        theta,
        σ,
        df,
        contig_lengths,
    )
end
m = exponential_density(df, contig_lengths)
chain = sample(m, NUTS(), 1000;)
DataFrame(summarize(chain))
