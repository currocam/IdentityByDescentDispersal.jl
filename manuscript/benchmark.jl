# This code reproduces the example provided in the manuscript
# and saves the results for reproducibility
using CSV, DataFrames, Turing, StatsBase, StatsPlots
using NCDatasets, ArviZ, IdentityByDescentDispersal
using Random, Plots, JLD2
Random.seed!(1000)
data = load("docs/data/constant_density.jld2")
# Pretend we have a CSV file
CSV.write("ibd_dispersal_data.csv", data["df"])
# Example code from the manuscript
df = CSV.read("ibd_dispersal_data.csv", DataFrame)
contig_lengths = [1.0]
@model function constant_density(df, contig_lengths)
    De ~ Uniform(0, 1e8)
    σ ~ Uniform(0, 1e8)
    cll = composite_loglikelihood_constant_density(De, σ, df, contig_lengths)
    Turing.@addlogprob! cll
end

benchmark_results = DataFrame(iteration = Int[], runtime_seconds = Float64[])
for i = 1:5
    start_time = time()
    # Multi-start optimization: run MLE from 5 different starting points and choose the best
    m = constant_density(df, contig_lengths)
    best_mle = nothing
    best_loglikelihood = -Inf

    for i = 1:5
        mle_candidate = maximum_likelihood(m; lb = [0.0, 0.0], ub = [1e8, 1e8])
        ll = mle_candidate.lp

        if ll > best_loglikelihood
            best_loglikelihood = ll
            best_mle = mle_candidate
        end
    end
    mle_estimate = best_mle
    coefs = mle_estimate |> coeftable |> DataFrame
    end_time = time()
    runtime = end_time - start_time
    push!(benchmark_results, (i, runtime))
    println("  Runtime: $(round(runtime, digits=3)) seconds")

end
# Save benchmark results to CSV
CSV.write("manuscript/mle_benchmark_results.csv", benchmark_results)
println("Mean runtime: $(round(mean(benchmark_results.runtime_seconds), digits=3)) seconds")
println("Std runtime: $(round(std(benchmark_results.runtime_seconds), digits=3)) seconds")
