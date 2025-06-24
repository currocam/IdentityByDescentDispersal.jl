```@meta
EditURL = "simulate_constant_density_ibd_detection.jl"
```

# Constant density simulation

````julia
using PyCall, Random, CSV, DataFrames, StatsBase, PrettyTables
using IdentityByDescentDispersal
````

## Forward-in-time simulation
We will use SLiM to simulate a constant density population living in a 2D torus and `tskit` and tree-sequence recording to analyze the ground truth of IBD blocks.

````julia
run(`slim --version`)
````

````
Process(`slim --version`, ProcessExited(0))
````

Recall that the average distance between the "mother" and father are not equal. The distance in one axis to the father is modeled as
$Z_{\text{father}} = N(0, SM) + N(0, SD)$
where $SM$ is the mate choice kernel and $SD$ is the dispersal rate of the offspring. The distance to the mother is modeled as
$Z_{\text{mother}} = N(0, SD)$

````julia
seed = 1000
D = 200 # Global local-density
SD = 0.1 # Dispersal rate of the offspring
SM = 0.01 # Mate choice kernel
outpath = "s$(seed).trees"
run(
    `slim -s $seed -d D=$D -d SD=$SD -d SM=$SM -d OUTPATH="\"$outpath\"" constant_density.slim`,
)
````

````
Process(`slim -s 1000 -d D=200 -d SD=0.1 -d SM=0.01 -d 'OUTPATH="s1000.trees"' constant_density.slim`, ProcessExited(0))
````

## Data preprocessing
We analyze the simulated tree-sequence using the python library `tskit`.

````julia
tskit = pyimport("tskit")
ts = tskit.load(outpath);
````

We perform recapitation to ensure samples are fully coalesced

````julia
pyslim = pyimport("pyslim")
rts = pyslim.recapitate(ts, ancestral_Ne=D, recombination_rate=1e-8);
````

````
/Users/currocampuzano/.julia/conda/3/aarch64/lib/python3.12/site-packages/pyslim/slim_metadata.py:874: UserWarning: This tree sequence is not the current SLiM format, so some operations may not work. Use `pyslim.update( )` to update the tree sequence.
  warnings.warn(
/Users/currocampuzano/.julia/conda/3/aarch64/lib/python3.12/site-packages/msprime/ancestry.py:1290: TimeUnitsMismatchWarning: The initial_state has time_units=ticks but time is measured in generations in msprime. This may lead to significant discrepancies between the timescales. If you wish to suppress this warning, you can use, e.g., warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
  sim = _parse_sim_ancestry(

````

We take a sample of 100 diploid individuals.

````julia
rng = MersenneTwister(seed)
n_samples = 100
sampled = randperm(rng, ts.num_individuals)[1:n_samples] .- 1
nodes = reduce(vcat, [ts.individual(i).nodes for i in sampled])
ts = ts.simplify(samples = nodes);
````

We overlay neutral mutations

````julia
msprime = pyimport("msprime")
mts = msprime.sim_mutations(ts, rate=1e-7, random_seed = seed);
````

Create VCF file

````julia
n_dip_indv = Int(ts.num_samples / 2)
indv_names = [ "tsk_$(i)indv" for i in 0:(n_dip_indv - 1) ]
outvcf = "s$(seed).vcf"
py"""
with open($outvcf, "w") as vcf_file:
    ts = $mts
    ts.write_vcf(vcf_file, individual_names=$indv_names)
"""
````

Run IBD detection software
wget https://faculty.washington.edu/browning/hap-ibd.jar

````julia
run(`curl -o hap-ibd.jar https://faculty.washington.edu/browning/hap-ibd.jar`)
````

````
Process(`curl -o hap-ibd.jar https://faculty.washington.edu/browning/hap-ibd.jar`, ProcessExited(0))
````

Create a dummy genetic map in Plink format

````julia
mapfile = "s$(seed).plink.map"
open(mapfile, "w") do io
    print(io, "1\trs\t0\t1\t\n1\trs\t100\t100000000\t\n")
end
````

Execute HapIBD software

````julia
run(`java -jar hap-ibd.jar gt=$outvcf map=$mapfile out=s$(seed)`)
````

````
Process(`java -jar hap-ibd.jar gt=s1000.vcf map=s1000.plink.map out=s1000`, ProcessExited(0))
````

Preprocess detected IBD blocks

````julia
run(`curl -o merge-ibd-segments.jar https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.17Jan20.102.jar`)
postprocessed_file = "s$(seed).postprocessed.ibd"
run(`bash -c "gunzip -c s$(seed).ibd.gz | java -jar merge-ibd-segments.jar $outvcf $mapfile 0.6 1 > $postprocessed_file"`)
````

````
Process(`bash -c 'gunzip -c s1000.ibd.gz | java -jar merge-ibd-segments.jar s1000.vcf s1000.plink.map 0.6 1 > s1000.postprocessed.ibd'`, ProcessExited(0))
````

Read IBD blocks
TO-DO: I'm not really sure what that the SCORE column is...

````julia
colnames = [
    "ID1", "HAP1", "ID2", "HAP2", "CHR", "START", "END", "SCORE", "LENGTH"
    ]
df_ibds = CSV.read(postprocessed_file, DataFrame; header = colnames);
df_ibds.span = df_ibds.LENGTH ./ 100;
````

Then, we compute pairwise distances across individuals (in the torus)

````julia
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
    df = DataFrame(ID1 = String[], ID2 = String[], distance = Float64[])
    for i = 1:n
        for j = (i+1):n
            push!(df, (indv_names[i], indv_names[j], dist_matrix[i, j]))  # 0-based ID
        end
    end
    df
end;
````

The package IdentityByDescentDispersal.jl has a function to preprocess both `DataFrames` into a single long DataFrame:
Here, we use the default bins used Ringbauer et. al. in their paper.

````julia
bins, min_threshold = default_ibd_bins()
df_preprocessed = preprocess_dataset(df_ibds, df_dist, bins, min_threshold);
````

In order to speedup the inference we are going to also bin the distances:

````julia
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
````

## Inference
Recall that we compare the MLE estimate of the density not with the "global" local density, but with the local density.
This is because, under this regime, individuals are clumped together and therefore experiment a higher density (e.g there are not uniformly distributed).

````julia
local_density = ts.metadata["SLiM"]["user_metadata"]["D"][1]
dispersal_rate = ts.metadata["SLiM"]["user_metadata"]["SIGMA"][1]
ground_truth = local_density, dispersal_rate
````

````
(253.30778276266545, 0.07088723439378913)
````

Finally, we can compute the MLE estimate and compare it with the ground truth:

````julia
using Turing
@model function constant_density(df, contig_lengths)
    D ~ Uniform(0, 1000)
    σ ~ Uniform(0, 1)
    Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
end
contig_lengths = [1.0]
m = constant_density(df2, contig_lengths);
mle_estimate = maximum_likelihood(m)
coef_table = mle_estimate |> coeftable |> DataFrame
pretty_table(coef_table, backend = Val(:text))
````

````
┌────────┬───────────┬────────────┬─────────┬─────────────┬───────────┬───────────┐
│   Name │     Coef. │ Std. Error │       z │    Pr(>|z|) │ Lower 95% │ Upper 95% │
│ String │   Float64 │    Float64 │ Float64 │     Float64 │   Float64 │   Float64 │
├────────┼───────────┼────────────┼─────────┼─────────────┼───────────┼───────────┤
│      D │   698.378 │     54.905 │ 12.7198 │   4.593e-37 │   590.766 │   805.989 │
│      σ │ 0.0763682 │ 0.00375581 │ 20.3334 │ 6.51756e-92 │  0.069007 │ 0.0837295 │
└────────┴───────────┴────────────┴─────────┴─────────────┴───────────┴───────────┘

````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
