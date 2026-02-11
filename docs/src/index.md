```@meta
CurrentModule = IdentityByDescentDispersal
```

# IdentityByDescentDispersal.jl

`IdentityByDescentDispersal.jl` is a Julia package for estimating effective population densities and dispersal rates from observed spatial patterns of IBD shared blocks. It provides the building blocks for performing likelihood-based inference and integrates seamlessly with other statistical libraries in the Julia ecosystem, such as [Turing.jl](https://turinglang.org).

## Installation

You can install this package by running:

```julia
import Pkg
Pkg.add("IdentityByDescentDispersal")
```

## Documentation

This website contains the documentation for [IdentityByDescentDispersal](https://github.com/currocam/IdentityByDescentDispersal.jl).

The documentation is organised into several sections:

- [Theory overview](@ref): An informal description of the theory behind the package, which covers the main concepts (e.g. dispersal rates, IBD blocks) and assumptions of the model.
- [Basic usage](@ref): A comprehensive guide to the functions provided by the package, to be used both for inference and prediction.
- [Inference and model evaluation](@ref): A complete inference example with simulated data and a set of recommendations for automated data processing (using a `Snakemake` pipeline) and model evaluation.
- [API reference](@ref): A complete reference of the functions exported by the package.
