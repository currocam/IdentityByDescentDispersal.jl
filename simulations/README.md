In addition to tests that assert the correctness of the implementation, we also include the necessary code to perform continuous-in-space simulations using SLiM and tree-sequence recording. The goal of these simulations is to generate realistic datasets to test the method.

To this date, it is not clear in which regimes the approximation proposed by Ringbauer *et al.* breaks down. Similarly, it is unclear to what extent errors in the detection of IBD blocks from phased data affect inference. Researchers should use simulated data that matches their system as part of their research workflow to assess and possibly adapt their expectations.

To this end, we include two basic scripts that can be used as a template by others.

- [Compute MLE estimate of a constant population with the ground truth IBD blocks](simulate_constant_density_ground_truth.md)
- [Compute MLE estimate of a constant population with IBD blocks detected with HapIBD and post-processed with merge-ibd.jar](simulate_constant_density_ibd_detection.md)

## Habitat suitability

The assumption that individuals disperse according to a Brownian walk is usually no more than a convenient modelling choice introduced for mathematical tractability. It is of great interest to the practitioner to assess to what extent deviations from these modelling assumptions affect inference. As an example, we include the [simulate_habitat_suitability_ground_truth.jl](simulate_habitat_suitability_ground_truth.jl) file, which can be used or adapted to study the effect of habitat patches on the accuracy of the inference.
