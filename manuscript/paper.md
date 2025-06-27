---
title: "IdentityByDescentDispersal.jl: Inferring dispersal rates with identity-by-descent blocks"
tags:
  - Julia
  - genomics
authors:
  - name: Francisco Campuzano Jiménez
    orcid: 0000-0001-8285-9318
    equal-contrib: true
    affiliation: 1
affiliations:
  - name: University of Antwerp, Belgium
    index: 1
    ror: 008x57b05
date: today
bibliography: paper.bib
---

# Summary

The population density and per-generation dispersal distance of an organism are central parameters in the study of evolution and ecology. The dispersal rate is particularly interesting for conservation management of fragmented or invasive species [@driscoll_trajectory_2014]. There is a growing interest in developing statistical methods that exploit the increasingly available genetic data to estimate such quantities [@rousset_genetic_1997;@ringbauer_inferring_2017; @smith_dispersal_2023;@smith_dispersenn2_2023].

# Statement of need

@ringbauer_inferring_2017 proposed an inference scheme that estimates effective population density and effective dispersal rate from shared identity-by-descent blocks. Despite their promising results, there is to this date no general-purpose software implementation of their method.

In order to make the proposed inference scheme available to the broader audience of evolutionary biologists and conservation scientists, we present `IdentityByDescentDispersal.jl`, a Julia package with an efficient and easy-to-use implementation of the method. This package implements the core equations proposed by @ringbauer_inferring_2017 and can be used to perform maximum-likelihood estimation or Bayesian inference via a composite likelihood.

In previous work, inference was restricted to a limited family of functions for the population effective density in the form $D_e(t) = Dt^{-\beta}$ for which the theory was analytically tractable. This implementation makes two major software contributions. First, it takes advantage of the powerful Julia ecosystem and the implementation of @geoga_fitting_2022 to provide an implementation of the composite-likelihood that is 100% compatible with automatic differentiation (AD), even with respect to $\beta$. Second, we also provide an AD-compatible interface to compute composite-likelihoods for arbitrary functions of $D_e(t)$ that we evaluate numerically via Gaussian quadrature rules.

`IdentityByDescentDispersal.jl` is compatible with standard gradient-based parameter estimation software, which is typically more efficient and yields better convergence. It also enables the fitting of a broader range of population-density functions that might be more appropriate to capture the characteristics of the system. Because of these advantages and its intuitive interface, we believe it will encourage a broader audience to adopt the inference scheme proposed by @ringbauer_inferring_2017.

# Overview

`IdentityByDescentDispersal.jl` contains two main types of functions. The first type of functions shares the prefix `expected_ibd_blocks`. They allow you to calculate the expected density of identity-by-descent blocks per pair of individuals and per unit of block length for various demographic models by solving (numerically or analytically) equation \autoref{eq:1} [@ringbauer_inferring_2017].

$$
\mathbb{E}[N_L | r, \theta] = \int_0^\infty G4t^2 \exp(-2Lt) \cdot \Phi(t | r, \theta) \,dt
\label{eq:1}
$$

where $G$ is the length of the genome (Morgans), $t$ is time (generations in the past), $L$ is the length of the block (Morgans), and $\Phi(t| r, \theta)$ is the instantaneous coalesce rate at time $t$ of two homologous loci that are separated $r$ units apart in the present (at time $t=0$) under the demography specified by $\theta$. A slightly more complicated expression that accounts for chromosomal edges and diploidy is the default in `IdentityByDescentDispersal.jl`.

The second type of function shares the prefix `composite_loglikelihood`. They directly compute the composite likelihood of the data by assuming the observed number of shared identity-by-descent blocks shared by a pair of individuals that are $r$ units apart whose length fall in a small bin $[L, L+\Delta L]$ follows a Poisson distribution with mean $\lambda = \mathbb{E}[N_L | r, \theta] \Delta L$.

`IdentityByDescentDispersal.jl` allows for three different parameterizations of effective population density: a constant density, a power-density, and a user-defined density (see Table \autoref{tab:tab1}).

| Function suffix    | $D_e(t)$ formula     | Parameters                | Solver       |
| ------------------ | -------------------- | ------------------------- | ------------ |
| `constant_density` | $D_e(t)=D$           | $D,\ \sigma$              | Analytically |
| `power_density`    | $D_e(t)=Dt^{-\beta}$ | $D,\ \beta,\ \sigma$      | Analytically |
| `custom`           | User-defined         | User-defined and $\sigma$ | Numerically  |

: `IdentityByDescentDispersal.jl` functions support three different parameterizations that are indicated by their respective suffix. \label{tab:tab1}

First, we provide a simulation template in SLiM for forward-in-time population genetics simulation in a continuous space with tree-sequence recording [@haller_slim_2023; @haller_tree-sequence_2019]. This template is available at XXX and might be used to asses model assumptions, guide empirical analysis, and perform simulation-based calibration. Assessing the performance of the data with synthetic datasets is a crucial step as it is known that errors in the detection of identity-by-descent blocks are common [@browning_identity_2012].

Second, we have also implemented a bioinformatics pipeline that carries out a complete analysis from detecting identity-by-descent blocks to finding the maximum-likelihood estimate of the effective population density and the effective dispersal rate. It is shared as a Snakemake pipeline, a popular bioinformatics workflow management tool [@molder_sustainable_2021]. It detects identity-by-descent blocks using HapIBD [@zhou_fast_2020], post-processes them with Refined IBD [@browning_improving_2013], and produces a CSV file compatible for subsequent analysis with `IdentityByDescentDispersal.jl` via the `preprocess_dataset` function.

# Citations

# Acknowledgements

# References
