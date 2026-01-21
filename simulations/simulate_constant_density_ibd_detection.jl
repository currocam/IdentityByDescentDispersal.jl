# # Constant density simulation
using PyCall, Random, CSV, DataFrames, StatsBase, PrettyTables
using IdentityByDescentDispersal
# ## Forward-in-time simulation
run(`slim --version`)

# Recall that the average distance between the "mother" and father are not equal. The distance in one axis to the father is modeled as
# $Z_{\text{father}} = N(0, SM) + N(0, SD)$
# where $SM$ is the mate choice kernel and $SD$ is the dispersal rate of the offspring. The distance to the mother is modeled as
# $Z_{\text{mother}} = N(0, SD)$
seed = 1000
NE = 200 # Number of individuals
SD = 0.1 # Dispersal rate of the offspring
SM = 0.01 # Mate choice kernel
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
outpath = "s$(seed).trees"
run(
    `slim -s $seed -d NE=$NE -d SD=$SD -d SM=$SM -d OUTPATH="\"$outpath\"" constant_density.slim`,
)

# ## Data preprocessing
# We analyze the simulated tree-sequence using the python library `tskit`.
tskit = pyimport("tskit")
ts = tskit.load(outpath);

# We perform recapitation to ensure samples are fully coalesced
pyslim = pyimport("pyslim")
rts = pyslim.recapitate(ts, ancestral_Ne = NE, recombination_rate = 1e-8);

# We take a sample of 100 diploid individuals.
rng = MersenneTwister(seed)
n_samples = 100
sampled = randperm(rng, ts.num_individuals)[1:n_samples] .- 1
nodes = reduce(vcat, [ts.individual(i).nodes for i in sampled])
ts = ts.simplify(samples = nodes);

# We overlay neutral mutations
msprime = pyimport("msprime")
mts = msprime.sim_mutations(ts, rate = 1e-7, random_seed = seed);

# Create VCF file
n_dip_indv = Int(ts.num_samples / 2)
indv_names = ["tsk_$(i)indv" for i = 0:(n_dip_indv-1)]
outvcf = "s$(seed).vcf"
py"""
with open($outvcf, "w") as vcf_file:
    ts = $mts
    ts.write_vcf(vcf_file, individual_names=$indv_names)
"""

# Run IBD detection software
# wget https://faculty.washington.edu/browning/hap-ibd.jar
run(`curl -o hap-ibd.jar https://faculty.washington.edu/browning/hap-ibd.jar`)
# Create a dummy genetic map in Plink format for recombination rate 1e-8
mapfile = "s$(seed).plink.map"
open(mapfile, "w") do io
    print(io, "1\trs\t0\t1\t\n1\trs\t100\t100000000\t\n")
end

# Execute HapIBD software
run(
    `java -Duser.language=en -Duser.country=US -jar hap-ibd.jar gt=$outvcf map=$mapfile out=s$(seed)`,
)

# Preprocess detected IBD blocks
run(
    `curl -o merge-ibd-segments.jar https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar`,
)
postprocessed_file = "s$(seed).postprocessed.ibd"
run(
    `bash -c "gunzip -c s$(seed).ibd.gz | java -Duser.language=en -Duser.country=US -jar merge-ibd-segments.jar $outvcf $mapfile 0.6 1 > $postprocessed_file"`,
)

# Read IBD blocks
# TO-DO: I'm not really sure what that the SCORE column is...
colnames = ["ID1", "HAP1", "ID2", "HAP2", "CHR", "START", "END", "SCORE", "LENGTH"]
df_ibds = CSV.read(postprocessed_file, DataFrame; delim = "\t", header = colnames);
df_ibds.span = df_ibds.LENGTH ./ 100;

# Then, we compute pairwise distances across individuals
df_dist = let
    points = reduce(hcat, [collect(row) for row in ts.individual_locations])'
    dist_matrix = euclidean_distance(points)
    n = size(points, 1)
    df = DataFrame(ID1 = String[], ID2 = String[], distance = Float64[])
    for i = 1:n
        for j = (i+1):n
            push!(df, (indv_names[i], indv_names[j], dist_matrix[i, j]))  # 0-based ID
        end
    end
    df
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
# Recall that we compare the MLE estimate of the density not with the "global" local density, but with the local density.
# This is because, under this regime, individuals are clumped together and therefore experience a higher density (e.g., they are not uniformly distributed).
local_density = ts.metadata["SLiM"]["user_metadata"]["D"][1]
dispersal_rate = ts.metadata["SLiM"]["user_metadata"]["SIGMA"][1]
ground_truth = local_density, dispersal_rate

# Finally, we can compute the MLE estimate and compare it with the ground truth:
using Turing
@model function constant_density(df, contig_lengths)
    D ~ Uniform(0, 1000)
    σ ~ Uniform(0, 100)
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
mle_estimate = maximum_likelihood(m)
coef_table = mle_estimate |> coeftable |> DataFrame
pretty_table(coef_table)
