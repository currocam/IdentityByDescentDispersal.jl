# IdentityByDescentDispersal

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://currocam.github.io/IdentityByDescentDispersal.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://currocam.github.io/IdentityByDescentDispersal.jl/dev)
[![Test workflow status](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Docs.yml?query=branch%3Amain)

# Getting started

This package provides an efficient implementation of the inference scheme proposed by [H. Ringbauer, G. Coop and N. H. Barton (2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340342/) to estimate the mean dispersal rate and the effective population density of a population.

The package is implemented in the [Julia programming language](https://julialang.org/) and designed to be used from within a julia session. It integrates seamlessly with other statistical libraries in the julia ecosystem such as [Turing.jl](https://turinglang.org). However, we also provide an automated [Snakemake pipeline](https://snakemake.github.io) to perform a complete analysis: from detecting and post-processing IBD blocks to finding a preliminary maximum likelihood estimate.

You can install this package by running:
```julia
import Pkg
Pkg.add("IdentityByDescentDispersal")
```

This package provides the building blocks for performing likelihood-based inference of effective population densities and mean effective dispersal rates. Together with [Turing.jl](https://turinglang.org), it offers a flexible interface for fitting Bayesian or maximum-likelihood models.

```julia
using IdentityByDescentDispersal
using CSV, DataFrames, Turing, StatsPlots
df = CSV.read("ibd_dispersal_data.csv", DataFrame)
contig_lengths = [1.0] # in Morgans
@model function constant_density(df, contig_lengths)
  D ~ LogNormal(1, 1) # Effective population density
  σ ~ Exponential(1) # Mean effective dispersal rate
  # ⬇️ Composite log-likelihood function provided by IdentityByDescentDispersal.jl
  Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
end
chain = sample(m, NUTS(), 1000)
```

Please refer to the [documentation](https://currocam.github.io/IdentityByDescentDispersal.jl) of the package for a detailed description of its functionality and recommended usage.

Simulation of synthetic datasets plays a major role in statistical inference and model validation. In addition to the inference machinery, this package also provides a set of recipes for developing advanced forward-in-time population genetics simulations in a continuous space. More information can be found in the [simulations subdirectory](simulations/README.md).


## Community Guidelines

`IdentityByDescentDispersal.jl` is an open-source project and contributions are welcome. Users are encouraged to report bugs, request features, or ask questions by opening a GitHub issue.
