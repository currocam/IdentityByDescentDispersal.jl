```@meta
CurrentModule = IdentityByDescentDispersal
```

# IdentityByDescentDispersal

Documentation for [IdentityByDescentDispersal](https://github.com/currocam/IdentityByDescentDispersal.jl).

# Topic 1

```@example
using IdentityByDescentDispersal
using Plots

L = 0.01         # Block length threshold (in Morgans)
G = 0.01         # Genome length (in Morgans)
D = 1.0          # Effective population density
σ = 1.0          # Dispersal rate
r = range(0.01, 25.0, length=200)  # Distances

plot(
    r,
    expected_ibd_blocks_constant_density.(r, 1.0, 0.5, L, G),
    xlabel = "Distance (r)",
    ylabel = "E[IBD blocks]",
    label="D=1.0, σ=0.5"
)
plot!(
    r,
    expected_ibd_blocks_constant_density.(r, 2.0, 0.5, L, G),
    label="D=2.0, σ=0.5"
)
plot!(
    r,
    expected_ibd_blocks_constant_density.(r, 1.0, 0.8, L, G),
    label="D=1.0, σ=0.8"
)
```

Additionally, we suppport user-define trajectories of effective density.

```@example
using IdentityByDescentDispersal
using Plots

# Parameters
L = 0.01         # Block length threshold (in Morgans)
G = 0.01         # Genome length (in Morgans)
σ = 0.5          # Dispersal rate
θ = [1.0, 0.5, 2π]  # Parameters for De(t): D₀, a, ω
r_values = range(0.01, 25.0, length=200)  # Distances
t_values = range(0.0, 10.0, length=200)   # Time for plotting D_e(t)

# User-defined effective density function
function De(t, θ)
    D₀, a, ω = θ
    D₀ * (1 + a * sin(ω * t))
end

# Compute expected IBD blocks across r values
ibd_values = [expected_ibd_blocks_custom(r, De, θ, σ, L, G) for r in r_values]
density_values = [De(t, θ) for t in t_values]

# Plotting
plot(
    plot(t_values, density_values,
        xlabel="Time (t)", ylabel="Effective Density Dₑ(t)",
        label="θ = [D₀=1, a=0.5, ω=2π]",
        title="Effective Density Trajectory",
    ),
    plot(r_values, ibd_values,
        xlabel="Distance (r)", ylabel="E[IBD blocks]",
        label="σ = 0.5, L = 0.01, G = 0.01",
        title="Expected Number of IBD Blocks",
    ),
    layout = (2, 1), size=(600, 800)
)
```
