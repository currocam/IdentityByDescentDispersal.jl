# Simulations

To this date, it is not clear in which regimes the approximation proposed by Ringbauer breaks, nor the expected amount of error introduced when detecting IBD blocks from phased data. We advise researchers to simulate data as part of their research workflow.

To facilitate this, we provide a set of simulations of population genetics in a continuous space that use forward-in-time simulations (with SLiM), tree-sequence recording, and coalescent simulation (Haller et al., 2019; Baumdicker et al., 2022; Haller & Messer, 2023).

## Simulation-based calibration

Simulation-based calibration (using power-composite likelihoods) is an unexplored and interesting approach. Feel free to contact me if you are interested in this topic.

## Basic recipes

To this end, we include two basic scripts that can be used as a template by others.

- [Compute MLE estimate of a constant population with the ground truth IBD blocks](simulate_constant_density_ground_truth.md)
- [Compute MLE estimate of a constant population with IBD blocks detected with HapIBD and post-processed with merge-ibd.jar](simulate_constant_density_ibd_detection.md)

## Advanced recipes (WIP)

Additionally, we provide an (incomplete) set of more advanced recipes.

### Habitat suitability

The assumption that individuals disperse according to a Brownian walk is usually no more than a convenient modelling choice introduced for mathematical tractability. It is of great interest to the practitioner to assess to what extent deviations from these modelling assumptions affect inference. As an example, we include the [simulate_habitat_suitability_ground_truth.jl](simulate_habitat_suitability_ground_truth.jl) file, which can be used or adapted to study the effect of habitat patches on the accuracy of the inference.
