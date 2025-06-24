```@meta
EditURL = "tutorial.jl"
```

# Getting started

This tutorial demonstrates basic usage of the package.

````@example tutorial
using IdentityByDescentDispersal
````

## Expected number of identity-by-descent blocks
The main feature of this package is to compute the expected number of identity-by-descent blocks under different demographic scenarios.
For example, we might explore how fast the number of blocks decay as a function of physical distance between pairs of individuals.

````@example tutorial
using Plots

L = 0.01         # Block length threshold (in Morgans)
G = 0.01         # Genome length (in Morgans)
D = 1.0          # Effective population density
σ = 1.0          # Dispersal rate
r_values = range(0.01, 25.0, length=200);  # Distances

plot(
    r_values,
    expected_ibd_blocks_constant_density.(r_values, 1.0, 0.5, L, G),
    xlabel = "Distance (r)",
    ylabel = "E[IBD blocks]",
    title = "Constant effective density scenario",
    label = "D=1.0, σ=0.5"
)
plot!(
    r_values,
    expected_ibd_blocks_constant_density.(r_values, 2.0, 0.5, L, G),
    label="D=2.0, σ=0.5"
)
plot!(
    r_values,
    expected_ibd_blocks_constant_density.(r_values, 1.0, 0.8, L, G),
    label="D=1.0, σ=0.8"
)
````

In addition to constant effective densities and power density functions, we also support arbitrary effective density functions provided by the user via numerical integration.
Next, we illustrate how to define a custom function by considering a scenario with an oscillating effective density function.

````@example tutorial
function De(t, θ)
    D₀, a, ω = θ
    D₀ * (1 + a * sin(ω * t))
end

θ = [1.0, 0.5, 2π]  # Parameters for De(t): D₀, a, ω
t_values = range(0.0, 10.0, length=200)   # Time for plotting D_e(t)
ibd_values = [expected_ibd_blocks_custom(r, De, θ, σ, L, G) for r in r_values]
density_values = [De(t, θ) for t in t_values]
plot(
    plot(t_values, density_values,
        xlabel="Time (t)", ylabel="Effective Density Dₑ(t)",
        label="θ = [D₀=1, a=0.5, ω=2π]",
        title="Effective Density Trajectory",
    ),
    plot(r_values, ibd_values,
        xlabel="Distance (r)", ylabel="E[IBD blocks]",
        label="σ = 0.5, L = 0.01, G = 0.01",
        title="Expected Number of IBD Blocks",
    ),
    layout = (2, 1), size=(600, 800),
    title = "User-defined effective density function"
)
````

## Inference

Ringbauer proposed to do inference by assuming that the observed number of IBD blocks that a pair $r$ units apart that fall within a small bin $[L_i, L_i + \Delta L]$ follows a Poisson distribution with mean $E[N_{L_i}(r, \theta)] \Delta L$.

Therefore, the log probability of observing `Y` identity-by-descent blocks whose length fall in the bin $[L_i, L_i + \Delta L]$ from a pair of individuals that are $r$ units apart can be calculated simply as:

````@example tutorial
using Distributions
λ = expected_ibd_blocks_constant_density(
    0.2, # Distance,
    0.5, # Dispersal rate
    2.0, # Effective density
    0.01, # Block length threshold (in Morgans)
    1.0 # Genome length (in Morgans)
)
Y = 5 # Observed number of IBD blocks
logpdf(Poisson(λ), Y)
````

The inference scheme is based on a composite-likelihood that treats each bin from different pairs of individuals as independent. For computational reasons, it is also faster to aggregate observations that correspond to the same distance and bin.

The input data for this sort of analysis is:

- A `DataFrame` containing the length of identity-by-descent blocks shared across different individuals
- A `DataFrame` containing distances across pairs of individuals.
- A list of contig lengths to properly take into account chromosomal edges (in Morgans).

Preprocessing the data requires combining the different sources of information and binning the identity-by-descent blocks.
For this purpose, we provide a helper function that returns a `DataFrame` in a "long" format.

````@example tutorial
using DataFrames
ibd_blocks = DataFrame(
    ID1 = ["A", "B", "A", "C", "B"],
    ID2 = ["B", "A", "C", "A", "C"],
    span = [0.005, 0.012, 0.21, 0.10, 0.08] # Morgans
)
individual_distances = DataFrame(
    ID1 = ["A", "A", "B"],
    ID2 = ["B", "C", "C"],
    distance = [10.0, 20.0, 10.0] # (e.g., in kilometers)
)
contig_lengths = [1.0, 1.5] # Contig lengths (in Morgans)
bins = [0.01, 0.02, 0.03] # Right edges of bins (in Morgans)
min_length = 0.005 # Minimum IBD length to consider
df = preprocess_dataset(ibd_blocks, individual_distances, bins, min_length)
````

We can now use the df `DataFrame`to compute composite loglikelihoods of different parameters using the family of `composite_loglikelihood_*` functions.
The second major feature of this package is that functions are all compatible with automatic differentiation. Analytical solutions depend on a modified Bessel function of the second kind, which is not widely implemented in a way that allows for automatic differentiation .
Here, we rely on the `BesselK.jl` package for doing this. In addition, we rely on the system of `QuadGK.jl` for numerical integration compatible with automatic differentiation. This allows us to perform gradient-based optimization and interact with existing software.

We recommend using this package together with the `Turing.jl`, as it provides a convenient interface for Bayesian and frequentist inference.

### Maximum likelihood estimation

For example, we can estimate the parameters of a constant density model using maximum likelihood estimation:

````@example tutorial
using Turing, StatsBase, StatsPlots
@model function constant_density(df, contig_lengths)
    D ~ Uniform(0, 100)
    σ ~ Uniform(0, 20)
    Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
end
````

Generate a MLE estimate.

````@example tutorial
mle_estimate = maximum_likelihood(constant_density(df, contig_lengths))
DataFrame(coeftable(mle_estimate)) # computed from the Fisher information matrix
````

### Bayesian inference

Alternatively, we can do a standard Bayesian inference with any of the available inference algorithms such as NUTS. Here, we fit a power-density model using NUTS:

````@example tutorial
@model function power_density(df, contig_lengths)
    D ~ Truncated(Normal(100, 20), 0, Inf)
    σ ~ Truncated(Normal(1, 0.1), 0, Inf)
    β ~ Normal(0, 0.5)
    Turing.@addlogprob! composite_loglikelihood_power_density(D, β, σ, df, contig_lengths)
end
m = power_density(df, contig_lengths)
````

Sample 4 chains in a serial fashion.

````@example tutorial
chains = sample(m, NUTS(), MCMCSerial(), 1000, 4)
plot(chains)
````

More complex models can be built with cautious consideration of identifiability. For example, we can fit a piecewise exponential density model using a custom function and numerical integration.

````@example tutorial
function piecewise_D(t, parameters)
    De1, De2, alpha, t0 = parameters
    t <= t0 ? De1*exp(-alpha*t) : De2
end
@model function exponential_density(df, contig_lengths)
    De1 ~ Truncated(Normal(1000, 100), 0, Inf)
    De2 ~ Truncated(Normal(1000, 100), 0, Inf)
    α ~ Normal(0, 0.05)
    t0 ~ Truncated(Normal(500, 100), 0, Inf)
    σ ~ Truncated(Normal(1, 0.1), 0, Inf)
    theta = [De1, De2, α, t0]
    Turing.@addlogprob! composite_loglikelihood_custom(piecewise_D, theta, σ, df, contig_lengths)
end
m = exponential_density(df, contig_lengths)
chain = sample(m, NUTS(), 1000;)
DataFrame(summarize(chain))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
