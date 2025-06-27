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

In previous work, inference was restricted to a limited family of functions for the population effective density in the form $D_e(t) = Dt^{-\beta}$ for which the theory was analytically tractable. This implementation makes two major software contributions. First, it takes advantage of the powerful Julia ecosystem and the implementation of Geoga et al. (2022) to provide an implementation of the composite-likelihood that is 100% compatible with automatic differentiation (AD), even with respect to $\beta$. Second, we also provide an AD-compatible interface to compute composite-likelihoods for arbitrary functions of $D_e(t)$ that we evaluate numerically via Gaussian quadrature rules.

`IdentityByDescentDispersal.jl` is compatible with standard gradient-based parameter estimation software, which is typically more efficient and yields better convergence. It also enables the fitting of a broader range of population-density functions that might be more appropriate to capture the characteristics of the system. Because of these advantages and its intuitive interface, we believe it will encourage a broader audience to adopt the inference scheme proposed by @ringbauer_inferring_2017.

# Overview

# Citations

# Acknowledgements

# References
