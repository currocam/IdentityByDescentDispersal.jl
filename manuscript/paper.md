---
title: "IdentityByDescentDispersal.jl: Inferring dispersal rates with identity-by-descent blocks"
tags:
  - Julia
  - population-genetics
  - IBD
  - dispersal
  - spatial-genetics
  - coalescent-theory
authors:
  - name: Francisco Campuzano Jiménez
    orcid: 0000-0001-8285-9318
    affiliation: 1
    corresponding: true
  - name: Arthur Zwaenepoel
    orcid: 0000-0003-1085-2912
    affiliation: 1
  - name: Els Lea R De Keyzer
    orcid: 0000-0003-0924-0118
    affiliation: 1
  - name: Hannes Svardal
    orcid: 0000-0001-7866-7313
    affiliation: "1, 2"
affiliations:
  - name: University of Antwerp, Belgium
    index: 1
    ror: 008x57b05
  - name: Naturalis Biodiversity Center, Leiden, Netherlands
    index: 2
    ror: 0566bfb96
date: 8 July 2025
bibliography: paper.bib
---

# Summary

The population density and per-generation dispersal rate of a population are central parameters in the study of evolution and ecology. The dispersal rate is particularly relevant for conservation management of fragmented or invasive species [@driscoll_trajectory_2014]. There is a growing interest in developing statistical methods that exploit the increasingly available genetic data to estimate the effective population density and effective dispersal rate [@rousset_genetic_1997;@ringbauer_inferring_2017;@smith_dispersal_2023;@smith_dispersenn2_2023].

The distribution of recent coalescent events between individuals in space can be used to estimate such quantities through the distribution of identity-by-descent (IBD) blocks [@ringbauer_inferring_2017]. An IBD block is defined as a segment of DNA that has been inherited by a pair of individuals from a common ancestor without being broken by recombination. Here we present `IdentityByDescentDispersal.jl`, a Julia package for estimating effective population densities and dispersal rates from observed spatial patterns of IBD shared blocks.

# Statement of need

@ringbauer_inferring_2017 proposed an inference scheme for the estimation of effective population density and effective dispersal rate from shared IBD blocks. Despite their promising results, there is to this date no general-purpose software implementation of their method.

In order to make the inference approach available to the broader audience of evolutionary biologists and conservation scientists, we present `IdentityByDescentDispersal.jl`, a Julia [@bezanson_julia_2017] package with an efficient and easy-to-use implementation of the method. The package implements the core equations proposed by @ringbauer_inferring_2017 and can be used to perform composite likelihood-based inference using either maximum-likelihood estimation (MLE) or Bayesian inference.

The method of @ringbauer_inferring_2017 was limited to a family of functions for the change in effective population density over time of the form $D_e(t) = Dt^{-\beta}$, for which the theory was analytically tractable. In addition, in the paper describing the original approach, the authors used gradient-free optimization to calculate maximum likelihood estimates (MLEs). Our implementation makes two major software contributions. First, we admit composite likelihood calculations for arbitrary functions $D_e(t)$ by evaluating the relevant integrals numerically through Gaussian quadrature rules [@quadgk]. Second, our implementation takes advantage of the powerful Julia ecosystem and the work of @geoga_fitting_2022 to provide a version of the composite likelihood that is fully compatible with automatic differentiation (AD), including AD with respect to $\beta$. By having a fully AD-compatible composite likelihood, `IdentityByDescentDispersal.jl` can be used together with standard gradient-based optimization and sampling methods available in the Julia ecosystem, which are typically more efficient than gradient-free methods.

Lastly, our package comes with a template to simulate synthetic datasets and a pipeline for end-to-end analysis from VCF files to final estimates. We believe it will encourage a broader audience to adopt the inference scheme proposed by @ringbauer_inferring_2017, motivate further developments and expand its applications.

# Overview

`IdentityByDescentDispersal.jl` contains two main sets of functions. The first set has the prefix `expected_ibd_blocks` and allows users to calculate the expected density of IBD blocks per pair of individuals and per unit of block length for various demographic models by solving \autoref{eq:1}.

\begin{align}\label{eq:1}
\mathbb{E}[N_L | r, \theta] = \int_0^\infty G4t^2 \exp(-2Lt) \cdot \Phi(t | r, \theta) \,dt
\end{align}

where $G$ is the length of the genome (in Morgan), $t$ is time (generations in the past), $L$ is the length of the block (Morgan) and $r$ is the geographical distance in the present (at time $t=0$) between the two individuals. $\Phi(t| r, \theta)$ is the instantaneous coalescence rate at time $t$ of two homologous loci that are initially $r$ units apart under the demographic model with parameters $\theta$. A slightly more complicated expression that accounts for chromosomal edges and diploidy is the default in `IdentityByDescentDispersal.jl`.

The second set of functions has the prefix `composite_loglikelihood` and allows users to directly compute the composite likelihood of the data by assuming the observed number of IBD blocks whose lengths fall in a small bin $[L, L+\Delta L]$ and are shared by a pair of individuals $r$ units apart follows a Poisson distribution with mean $\lambda = \mathbb{E}[N_L | r, \theta] \Delta L$.

`IdentityByDescentDispersal.jl` allows for three different parameterizations of the effective population density function: a constant density, a power-density, and a user-defined density (see \autoref{tab:tab1}).

| Function suffix    | $D_e(t)$ formula     | Parameters                | Solver       |
| ------------------ | -------------------- | ------------------------- | ------------ |
| `constant_density` | $D_e(t)=D$           | $D,\ \sigma$              | Analytically |
| `power_density`    | $D_e(t)=Dt^{-\beta}$ | $D,\ \beta,\ \sigma$      | Analytically |
| `custom`           | User-defined         | User-defined and $\sigma$ | Numerically  |

: `IdentityByDescentDispersal.jl` functions support three different parameterizations that are indicated by their respective suffixes. \label{tab:tab1}

The Julia package is accompanied by two additional resources. First, we provide a simulation template in SLiM for forward-in-time population genetics simulation in a continuous space with tree-sequence recording [@haller_slim_2023; @haller_tree-sequence_2019]. This template can be used to assess model assumptions, guide empirical analysis, and perform simulation-based calibration. Assessing the performance of the method with synthetic datasets is a crucial step, as it is known that errors in the detection of IBD blocks are common [@browning_identity_2012] and that inferences based on composite likelihood tend to be overconfident, underestimating posterior uncertainty and yielding too narrow confidence intervals.

Second, we have also implemented a bioinformatics pipeline that carries out a complete analysis from detecting IBD blocks to finding the MLE of the effective population density and the effective dispersal rate. It is shared as a Snakemake pipeline, a popular bioinformatics workflow management tool [@molder_sustainable_2021]. It takes as input a set of phased VCF files, their corresponding genetic maps and a CSV file containing pairwise geographical distances between individuals. The pipeline detects IBD blocks using HapIBD [@zhou_fast_2020], post-processes them with Refined IBD [@browning_improving_2013] and produces a CSV file compatible with subsequent analysis with `IdentityByDescentDispersal.jl` via the `preprocess_dataset` function.

Both the SLiM simulation template and the Snakemake pipeline can be found in the GitHub repository at [https://github.com/currocam/IdentityByDescentDispersal.jl](https://github.com/currocam/IdentityByDescentDispersal.jl).

# Example

In this section, we demonstrate how `IdentityByDescentDispersal.jl` can be used together with the popular `Turing.jl` framework [@ge_turing_2018] using a dataset we simulate in the documentation. We analyze error-free IBD blocks shared by 100 diploid individuals from a constant-density population with parameters $D_{\text{true}}\approx 250$ diploids/km² and $\sigma_{\text{true}}\approx 0.071$ km/generation.

`IdentityByDescentDispersal.jl` has extensive documentation that covers the underlying theory behind the method, how to effectively simulate synthetic datasets, various demographic models, and inference algorithms. We refer the reader to the documentation for more details, which can be found at [https://currocam.github.io/IdentityByDescentDispersal.jl/](https://currocam.github.io/IdentityByDescentDispersal.jl/dev/).

Thanks to `Turing.jl`, we can perform Bayesian inference with a wide range of popular Monte Carlo algorithms. \autoref{fig:example} shows the estimated pseudo-posterior obtained through doing inference with the composite likelihood.

![Estimated pseudo-posterior obtained by doing inference with the composite likelihood. Although the pseudo-posterior is not well calibrated, it concentrates near the true values ($\mathbb E[D | \text{data}]\approx 281$ and $\mathbb E[\sigma | \text{data}]\approx 0.068$, respectively).\label{fig:example}](figures/nuts_constant_density.svg){ width=100% }

\autoref{fig:example} was generated by the following snippet of Julia code, which reads the processed data CSV from the provided Snakemake pipeline.

```julia
using CSV, DataFrames, Turing, StatsPlots, IdentityByDescentDispersal
df = CSV.read("ibd_dispersal_data.csv", DataFrame)
contig_lengths = [1.0]
@model function constant_density(df, contig_lengths)
    De ~ Truncated(Normal(1000, 100), 0, Inf)
    σ ~ InverseGamma(1, 1)
    Turing.@addlogprob! composite_loglikelihood_constant_density(
      De, σ, df, contig_lengths
    )
end
m = constant_density(df, contig_lengths)
chains = sample(m, NUTS(), MCMCThreads(), 1000, 4)
plot(chains)
```

We can also easily compute the MLEs of the same demographic model,

```julia
mle_estimate = maximum_likelihood(
  m; lb=[0.0, 0.0], ub=[1e8, 1e8]
)
coeftable(mle_estimate)
```

which estimates $D_{\text{MLE}}\approx 282$ diploids/km² (95% CI: 260–303) and $\sigma_{\text{MLE}}\approx 0.068$ km/generation (95% CI: 0.065–0.071). The 95% confidence interval is computed from the Fisher information matrix.

# Availability

`IdentityByDescentDispersal.jl` is a registered Julia package available through the official General registry. Its source code is hosted on GitHub at [https://github.com/currocam/IdentityByDescentDispersal.jl](https://github.com/currocam/IdentityByDescentDispersal.jl).

# Acknowledgements

We acknowledge financial support from the Research Foundation - Flanders (FWO). This work was supported by FWO-G0A9B24N (F.C.J, H.S), FWO-1272625N (A.Z., H.S) and FWO-12A8423N (E.L.R.D.K., H.S).

# References
