# IdentityByDescentDispersal

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://currocam.github.io/IdentityByDescentDispersal.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://currocam.github.io/IdentityByDescentDispersal.jl/dev)
[![Test workflow status](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/currocam/IdentityByDescentDispersal.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/currocam/IdentityByDescentDispersal.jl)
[![Lint workflow Status](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/currocam/IdentityByDescentDispersal.jl/actions/workflows/Docs.yml?query=branch%3Amain)

# Getting started

You can install this package by running:

```julia
import Pkg
Pkg.add(url="https://github.com/currocam/IdentityByDescentDispersal.jl")
```

Please refer to the [documentation](https://currocam.github.io/IdentityByDescentDispersal.jl) of the package for a detailed description of this package.

## Overview

This Julia package implements a method proposed by [Ringbauer et al. (2017)](10.1534/genetics.116.196220) to estimate effective densities and effective dispersal rates from identity-by-descent blocks.

This method assumes a diffusion approximation to trace genetic ancestry and leverages coalescent theory to derive analytical formulas for the expected density of identity-by-descent blocks shared between a pair of individuals that are $r$ units apart according to a certain demographic model $\theta$.

```julia
using IdentityByDescentDispersal
r = 0.5         # Distance between pair of diploid individuals
L = 0.01        # Block length threshold (in Morgans)
G = 1.0         # Genome length (in Morgans)
D = 1.0         # Effective population density
σ = 1.0         # Dispersal rate
expected_ibd_blocks_constant_density(r, D, σ, L, G)
```

```
1605.2440136343164
```

In addition to scenarios that lead to analytical solutions, we support arbitrary demographic models that are solved via efficient numerical integration using Gaussian-quadrature rules with [QuadGK.jl](https://juliamath.github.io/QuadGK.jl/stable/).

```julia
# Exponential density (backwards in time)
De(t, θ) = θ[1] * exp(-θ[2] * t)
α = -0.01 # Growth rate
expected_ibd_blocks_custom(r, De, [D, α], σ, L, G)
```

```
719.9535621865995
```

[Ringbauer et al. (2017)](10.1534/genetics.116.196220) proposed an inference scheme based on a composite likelihood. This package implements an improved version of this scheme that relies on automatic differentiation and that integrates well with Julia's ecosystem.

Analytical solutions to the expected density of identity-by-descent blocks depend on the modified second-kind Bessel function $K_v(x)$. Here we rely on the implementation of [Geoga et al. (2022)](https://arxiv.org/pdf/2201.00090) to compute fast and accurate derivatives with respect to both $v$ and $x$. Therefore, this package is fully compatible with automatic differentiation using [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) when using both analytical solutions and numerical approximations. This leads to gradient-based optimization that is generally more efficient when finding maximum-likelihood estimates and allows us to use standard Bayesian inference software that uses Hamiltonian Monte Carlo algorithms.

```julia
using CSV, DataFrames, Turing
df = CSV.read("ibd_dispersal_data.csv", DataFrame)
contig_lengths = [1.0] # in Morgans
@model function constant_density(df, contig_lengths)
    D ~ LogNormal(log(10), 1.0)
    σ ~ truncated(Normal(0, 0.5), 0, 1)
    Turing.@addlogprob! composite_loglikelihood_constant_density(D, σ, df, contig_lengths)
end
chain = sample(m, NUTS(), 1000)
mean(chain[:D]), mean(chain[:σ])
```

```
(23.70889539436043, 0.47267472059392046)
```

Please refer to the [documentation](https://currocam.github.io/IdentityByDescentDispersal.jl/dev/) for the complete reference of functions and detailed overview of this package and its usage to perform maximum-likelihood or Bayesian inference.

## Simulations

To this date, it is not entirely clear under which regimes the diffusion approximation breaks. We advise researchers to use simulations to adequate their expectations and guide their analysis using simulated datasets, especially when choosing the identity-by-descent detection software.

We provide a set of simulations of population genetics simulations in a continuous space that use forward-in-time simulations (with SLiM), tree-sequence recording, and coalescent simulation (Haller et al., 2019; Baumdicker et al., 2022; Haller & Messer, 2023). These simulations can be found in the [simulations directory](simulations/README.md) and can be used as templates.

Simulation-based calibration (using power-composite likelihoods) is an unexplored and interesting approach. Feel free to contact us if you are interested.

## Full analysis workflow

Inferring the effective density and dispersal rate from empirical datasets requires running a identity-by-descent detection software across different phased VCF files, post-processing them, and aggregating them appropriately. The function `preprocess_dataset`is designed to simplify this process.

In addition, we include in this repository a [Snakemake pipeline](Snakefile) that can be used to perform a complete analysis, from detecting IBD blocks using [HapIBD](https://github.com/browning-lab/hap-ibd), post-processing them with [Refined IBD](https://faculty.washington.edu/browning/refined-ibd.html), producing a CSV directly compatible with this package, and, optionally, finding a preliminary maximum likelihood estimate of the effective density and effective dispersal rate.

Such a pipeline can be configured with a `config.yaml` file and executed by running:

```bash
snakemake -s Snakefile --configfile config.yaml --cores 8 --sdm conda
```

You may find an example of a [configuration file](.test-workflow/config.yaml) and compatible data at the [.test-workflow subdirectory](.test-workflow).

## References

Baumdicker, F., Bisschop, G., Goldstein, D., Gower, G., Ragsdale, A. P., Tsambos, G., Zhu, S., Eldon, B., Ellerman, E. C., Galloway, J. G., Gladstein, A. L., Gorjanc, G., Guo, B., Jeffery, B., Kretzschumar, W. W., Lohse, K., Matschiner, M., Nelson, D., Pope, N. S., … Kelleher, J. (2022). Efficient ancestry and mutation simulation with msprime 1.0. Genetics, 220(3), iyab229. https://doi.org/10.1093/genetics/iyab229
Geoga, C. J., Marin, O., Schanen, M., & Stein, M. L. (2022). Fitting Matérn Smoothness Parameters Using Automatic Differentiation (Version 3). arXiv. https://doi.org/10.48550/ARXIV.2201.00090
Haller, B. C., Galloway, J., Kelleher, J., Messer, P. W., & Ralph, P. L. (2019). Tree-sequence recording in SLiM opens new horizons for forward-time simulation of whole genomes. Molecular Ecology Resources, 19(2), 552–566. https://doi.org/10.1111/1755-0998.12968
Haller, B. C., & Messer, P. W. (2023). SLiM 4: Multispecies Eco-Evolutionary Modeling. The American Naturalist, 201(5), E127–E139. https://doi.org/10.1086/723601
Ringbauer, H., Coop, G., & Barton, N. H. (2017). Inferring Recent Demography from Isolation by Distance of Long Shared Sequence Blocks. Genetics, 205(3), 1335–1351. https://doi.org/10.1534/genetics.116.196220
