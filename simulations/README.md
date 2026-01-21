In addition to tests that assert the correctness of the implementation, we also include the necessary code to perform continuous-in-space simulations using SLiM and tree-sequence recording. Their goal is to simulate realistic datasets to test the method.

To this date, it is not clear in which regimes the approximation proposed by Ringbauer breaks, nor the expected amount of error introduced when detecting IBD blocks from phased data. Researchers should simulate data as part of their research workflow that matches their system to adapt their expectations.

We include two basic scripts that can be adapted by others to match their case.

- [Compute MLE estimate of a constant population with the ground truth IBD blocks](simulate_constant_density_ground_truth.md)
- [Compute MLE estimate of a constant population with IBD blocks detected with HapIBD and post-processed with merge-ibd.jar](simulate_constant_density_ibd_detection.md)

## Habitat suitability

In practice, the assumption that individuals disperse according to a Brownian walk is not a modelling choice, but rather a statistical requirement for tractability. It is of great interest to the practitioner to address how such deviations from the model assumptions affect the inference. The [simulate_habitat_suitability_ground_truth.jl](simulate_habitat_suitability_ground_truth.jl) file can be used or adapted to study the effect of habitat patches on the accuracy of the inference.
