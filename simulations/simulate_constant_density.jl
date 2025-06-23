# # Constant density simulation
using PyCall, Random, DataFrames, StatsBase, MarkdownTables
using IdentityByDescentDispersal
# ## Forward-in-time simulation
# We will use SLiM to simulate a constant density population living in a 2D torus and `tskit` and tree-sequence recording to analyze the ground truth of IBD blocks.

run(`slim --version`)

# Recall that the average distance between the "mother" and father are not equal. The distance in one axis to the father is modeled as
# $Z_{\text{father}} = N(0, SM) + N(0, SD)$
# where $SM$ is the mate choice kernel and $SD$ is the dispersal rate of the offspring. The distance to the mother is modeled as
# $Z_{\text{mother}} = N(0, SD)$
seed = 1000
D = 200 # Global local-density
SD = 0.1 # Dispersal rate of the offspring
SM = 0.01 # Mate choice kernel
outpath = "s$(seed).trees"
run(
    `slim -s $seed -d D=$D -d SD=$SD -d SM=$SM -d OUTPATH="\"$outpath\"" constant_density.slim`,
)

# ## Data preprocessing
# We analyze the simulated tree-sequence using the python library `tskit`. Because IBD blocks decay very quickly with time and
# we have run the simulation for many generations (5000), we don't need to perform `recapitation`.
tskit = pyimport("tskit")
ts = tskit.load(outpath);

# We take a sample of 200 diploid individuals.
rng = MersenneTwister(seed)
n_samples = 100
sampled = randperm(rng, ts.num_individuals)[1:n_samples] .- 1
nodes = reduce(vcat, [ts.individual(i).nodes for i in sampled])
ts = ts.simplify(samples = nodes);

# Then, we compute pairwise distances across individuals (in the torus)
function torus_distance(points::AbstractMatrix{<:Real})
    coords = points[:, 1:2]
    n = size(coords, 1)
    dist_matrix = zeros(n, n)

    for i = 1:n
        for j = 1:n
            dx = abs(coords[i, 1] - coords[j, 1])
            dy = abs(coords[i, 2] - coords[j, 2])
            dx = min(dx, 1 - dx)
            dy = min(dy, 1 - dy)
            dist_matrix[i, j] = sqrt(dx^2 + dy^2)
        end
    end
    return dist_matrix
end
df_dist = let
    points = ts.individual_locations
    dist_matrix = torus_distance(points)
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
            (block.right - block.left) / 1e8 for pair in ibds for block in ibds.get(pair)
        ]
        DataFrame(ID1 = i, ID2 = j, span = spans)
    end
    vcat([blocks(i, j) for i = 0:(n_samples-1) for j = (i+1):(n_samples-1)]...)
end;

# The package IdentityByDescentDispersal.jl has a function to preprocess both `DataFrames` into a single long DataFrame:
# Here, we use the default bins used Ringbauer et. al. in their paper.
bins, min_threshold = default_ibd_bins()
df_preprocessed = preprocess_dataset(df_ibds, df_dist, bins, min_threshold);

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

df2 = let
    edges_distance = collect(0:0.05:1.0)
    df_preprocessed.DISTANCE_INDEX = cut(df_preprocessed.DISTANCE, edges_distance)
    df2 = combine(
        groupby(df_preprocessed, [:DISTANCE_INDEX, :IBD_LEFT, :IBD_RIGHT]),
        :NR_PAIRS => sum => :NR_PAIRS,
        :COUNT => sum => :COUNT,
        :DISTANCE => mean => :DISTANCE,
    )
end;

# ## Inference
# Finally, we can compute the MLE estimate and compare it with thr ground truth:
# Recall that we compare the MLE estimate of the density not with the "global" local density, but with the local density.
# This is because, under this regime, individuals are clumped together and therefore experiment a higher density (e.g there are not uniformly distributed).
local_density = ts.metadata["SLiM"]["user_metadata"]["D"][1]
dispersal_rate = ts.metadata["SLiM"]["user_metadata"]["SIGMA"][1]
ground_truth = local_density, dispersal_rate

using Turing
@model function constant_density(df, contig_lengths)
    D ~ Uniform(0, 1000)
    σ ~ Uniform(0, 1)
    Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
end
contig_lengths = [1.0]
m = constant_density(df2, contig_lengths);
mle_estimate = maximum_likelihood(m)
mle_estimate |> coeftable |> DataFrame |> markdown_table()
