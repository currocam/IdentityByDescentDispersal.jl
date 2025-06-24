using Literate
Literate.markdown(
    "simulate_constant_density_ground_truth.jl",
    ".";
    name = "simulate_constant_density_ground_truth",
    execute = true,
    softscope = true,
)

Literate.markdown(
    "simulate_constant_density_ibd_detection.jl",
    ".";
    name = "simulate_constant_density_ibd_detection",
    execute = true,
    softscope = true,
)
