# # Habitat suitability simulation
using PyCall, Random, DataFrames, StatsBase, PrettyTables
using IdentityByDescentDispersal
using Turing
# ## Forward-in-time simulation

run(`slim --version`)

function euclidean_distance(points)
    coords = points[:, 1:2]
    n = size(coords, 1)
    dist_matrix = zeros(n, n)
    for i = 1:n
        for j = 1:n
            dx = coords[i, 1] - coords[j, 1]
            dy = coords[i, 2] - coords[j, 2]
            dist_matrix[i, j] = sqrt(dx^2 + dy^2)
        end
    end
    return dist_matrix
end
# In order to speedup the inference we are going to also bin the distances:
function cut(values::AbstractVector{<:Real}, edges::AbstractVector{<:Real})
    n = length(values)
    bins = Vector{Float64}(undef, n)
    for i = 1:n
        x = values[i]
        bin_idx = findfirst(j -> edges[j] <= x < edges[j+1], 1:(length(edges)-1))
        @assert !isnothing(bin_idx)
        bins[i] = edges[bin_idx]  # label bin by left edge
    end
    return bins
end

# Define the simulation loop
function simulation(NE, SD, SM, seed)
    outpath = "habitat_suitability.trees"
    run(
        `slim -s $seed -d NE=$NE -d SD=$SD -d SM=$SM -d OUTPATH="\"$outpath\"" habitat_suitability.slim`,
    )
    # ## Data preprocessing
    # We analyze the simulated tree-sequence using the python library `tskit`. Because IBD blocks decay very quickly with time and
    # we have run the simulation for many generations (5000), we don't need to perform `recapitation`.
    tskit = pyimport("tskit")
    ts = tskit.load(outpath);

    # We take a sample of 200 diploid individuals.
    rng = MersenneTwister(seed)
    n_samples = 200
    sampled = randperm(rng, ts.num_individuals)[1:n_samples] .- 1
    nodes = reduce(vcat, [ts.individual(i).nodes for i in sampled])
    ts = ts.simplify(samples = nodes);

    df_dist = let
        points = reduce(hcat, [collect(row) for row in ts.individual_locations])'
        dist_matrix = euclidean_distance(points)
        n = size(points, 1)
        df = DataFrame(ID1 = Int[], ID2 = Int[], distance = Float64[])
        for i = 1:n
            for j = (i+1):n
                push!(df, (i - 1, j - 1, dist_matrix[i, j]))  # 0-based ID
            end
        end
        df
    end;

    # And the IBD blocks:
    df_ibds = let
        function blocks(i, j)
            ibds = ts.ibd_segments(
                between = [ts.individual(i).nodes, ts.individual(j).nodes],
                min_span = 1 / 100 * 4 * 1e8,
                store_pairs = true,
                store_segments = true,
            )
            spans = [
                (block.right - block.left) / 1e8 for pair in ibds for
                block in ibds.get(pair)
            ]
            DataFrame(ID1 = i, ID2 = j, span = spans)
        end
        vcat([blocks(i, j) for i = 0:(n_samples-1) for j = (i+1):(n_samples-1)]...)
    end;

    # The package IdentityByDescentDispersal.jl has a function to preprocess both `DataFrames` into a single long DataFrame:
    # Here, we use the default bins used by Ringbauer et. al. in their paper.
    bins, min_threshold = default_ibd_bins()
    df_preprocessed = preprocess_dataset(df_ibds, df_dist, bins, min_threshold);
    df2 = let
        h = fit(Histogram, df_preprocessed.DISTANCE, nbins = 20)
        df_preprocessed.DISTANCE_INDEX =
            cut(df_preprocessed.DISTANCE, collect(h.edges[1]))
        df2 = combine(
            groupby(df_preprocessed, [:DISTANCE_INDEX, :IBD_LEFT, :IBD_RIGHT]),
            :NR_PAIRS => sum => :NR_PAIRS,
            :COUNT => sum => :COUNT,
            :DISTANCE => mean => :DISTANCE,
        )
    end;
    # ## Inference
    # Recall that we compare the MLE estimate of the density not with the "global" local density, but with the local density.
    # This is because, under this regime, individuals are clumped together and therefore experience a higher density (e.g, they are not uniformly distributed).
    local_density = ts.metadata["SLiM"]["user_metadata"]["D"][1]
    dispersal_rate = ts.metadata["SLiM"]["user_metadata"]["SIGMA"][1]
    ground_truth = local_density, dispersal_rate
    # Finally, we can compute the MLE estimate and compare it with the ground truth:
    @model function constant_density(df, contig_lengths)
        D ~ Uniform(0, 500)
        σ ~ Uniform(0, 20.0)
        try
            Turing.@addlogprob! composite_loglikelihood_constant_density(
                D,
                σ,
                df,
                contig_lengths,
            )
        catch e
            @warn "Error in constant_density model: $e"
            Turing.@addlogprob! -Inf
        end
    end
    contig_lengths = [1.0]
    m = constant_density(df2, contig_lengths);
    mle_estimate = maximum_a_posteriori(m)
    DataFrame(
        Dict(
            :D_truth => ground_truth[1],
            :D_estimate => mle_estimate.values[:D],
            :σ_truth => ground_truth[2],
            :σ_estimate => mle_estimate.values[:σ],
        ),
    )
end

Ne = 300
seed = 1000
SM = 0.01
kernels_sd = collect(range(0.1, 3.0, 50))
using ProgressBars
simulations = [simulation(Ne, SD, SM, seed) for SD in ProgressBar(kernels_sd)]
df_sims = reduce(vcat, simulations)
using CSV
CSV.write("habitat_sims.csv", df_sims)


### Quick visualization of the results
using Plots
# Create scatter plot for σ
p1 = scatter(
    df_sims.σ_truth,
    df_sims.σ_estimate,
    label = "",
    xlabel = "True dispersal rate",
    ylabel = "MLE dispersal rate",
)
cor_sigma = cor(df_sims.σ_truth, df_sims.σ_estimate)
plot!(p1, x -> x, df_sims.σ_truth[1], df_sims.σ_truth[end], label = "")
title!(p1, "correlation(true σ, MLE σ): $(round(cor_sigma, digits=2))")

# How σ affects the estimated density
p2 = scatter(
    df_sims.σ_estimate,
    df_sims.D_estimate,
    label = "Maximum likelihood",
    xlabel = "Dispersal rate (σ)",
    ylabel = "Effective density (De)",
)
scatter!(df_sims.σ_truth, df_sims.D_truth, label = "Ground truth")

# Visualize the map used in the simulation
# Notice that SLiM defines the maps slightly different than how heatmap works
habitat_suitability = ones(100, 100)
habitat_suitability[36:65, 36:65] .= 0.05

# Visualize sampling points only for the last simulations
ts = pyimport("tskit").load("habitat_suitability.trees");
all_points = reduce(hcat, [collect(row) for row in ts.individual_locations])'
rng = MersenneTwister(seed)
n_samples = 200
sampled = randperm(rng, ts.num_individuals)[1:n_samples] .- 1
nodes = reduce(vcat, [ts.individual(i).nodes for i in sampled])
ts = ts.simplify(samples = nodes);
points = reduce(hcat, [collect(row) for row in ts.individual_locations])'
p3 = heatmap(
    habitat_suitability,
    aspect_ratio = :equal,
    color = :grays,
    axis = nothing,
    grid = false,
    colorbar_title = "Habitat suitability",
    xlim = (0, 100),
    ylim = (0, 100),
    title = "Simulation setup (De=300, 1x1 grid)",
    colorbar = true,
)
using KernelDensity
x = Float64.(all_points[:, 1]) * 100
y = Float64.(all_points[:, 2]) * 100
k = kde((x, y))
p4 = contourf(
    k.x / 100,
    k.y / 100,
    k.density / sum(k.density) * 300,
    aspect_ratio = :equal,
    colorbar = true,
    xlim = (0, 1),
    ylim = (0, 1),
    colorbar_title = "Density kernel",
    title = "Representative replicate (De=300, σ=0.34)",
)
scatter!(
    p4,
    points[1:end, 1],
    points[1:end, 2],
    label = "Sampled points (n=200)",
    color = :grey,
    markershape = :xcross,
)


plot(p1, p2, layout = (2, 1), size = (800, 600))

gr()
Plots.GR.beginprint("habitat_suitability_experiment.pdf")
plot(p3, size = (800, 600), dpi = 300)
plot(p4, size = (800, 600), dpi = 300)
plot(p1, size = (800, 600), dpi = 300)
plot(p2, size = (800, 600), dpi = 300)
Plots.GR.endprint()
