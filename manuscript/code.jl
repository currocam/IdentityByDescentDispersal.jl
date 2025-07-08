# This code reproduces the example provided in the manuscript
# and saves the results for reproducibility
using CSV, DataFrames, Turing, StatsBase, StatsPlots
using NCDatasets, ArviZ, IdentityByDescentDispersal
using Random, Plots, JLD2
Random.seed!(1234)
data = load("docs/data/constant_density.jld2")
# Pretend we have a CSV file
CSV.write("ibd_dispersal_data.csv", data["df"])
# Example code from the manuscript
df = CSV.read("ibd_dispersal_data.csv", DataFrame)
contig_lengths = [1.0]
@model function constant_density(df, contig_lengths)
    De ~ Truncated(Normal(500, 100), 0, Inf)
    σ ~ InverseGamma(1, 1)
    cll = composite_loglikelihood_constant_density(De, σ, df, contig_lengths)
    Turing.@addlogprob! cll
end
m = constant_density(df, contig_lengths)
mle_estimate = maximum_likelihood(m; lb = [0.0, 0.0], ub = [1e8, 1e8])
coefs = mle_estimate |> coeftable |> DataFrame
chains = sample(m, NUTS(), MCMCThreads(), 1000, 4)
expected_posterior = chains |> mean |> DataFrame
# Save the DataFrame to a CSV file
CSV.write("manuscript/data/expected_posterior_constant_density.csv", expected_posterior)
CSV.write("manuscript/data/mle_constant_density.csv", coefs)
# Save MCMC chains for reproducibility
idata = from_mcmcchains(chains)
to_netcdf(idata, "manuscript/data/nuts_constant_density.nc")
default(
    fontfamily = "Computer Modern",
    tickfontsize = 10,
    guidefontsize = 12,
    legendfontsize = 10,
    extra_plot_kwargs = :overwrite_legend_position => false,
)
plt = plot(chains; size = (1000, 700))
savefig(plt, "manuscript/figures/nuts_constant_density.svg")
savefig(plt, "manuscript/figures/nuts_constant_density.pdf")
savefig(plt, "manuscript/figures/nuts_constant_density.png")
