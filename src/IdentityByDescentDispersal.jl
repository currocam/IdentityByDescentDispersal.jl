module IdentityByDescentDispersal
using BesselK: adbesselk

"""
    x = safe_adbesselk(1, 1.0)
A wrapper around the `adbesselk` to compute the modified Bessel function of the second kind of degree that returns a NaN if computation fails.
"""
@inline function safe_adbesselk(v, x)
    try
        return adbesselk(v, x)
    catch
        return NaN
    end
end

"""
    ϕ(t::Real, r::Real, D_e::Function, σ::Real)
Computes the probability ``\\phi(t)`` that two homologous loci coalesce ``t`` generations ago according to

```math
\\phi(t) = \\frac{1}{2D_e(t)} \\frac{1}{4\\pi t \\sigma^2} \\exp(\\frac{-r^2}{4 t \\sigma^2})
```

Ringbauer, H., Coop, G. & Barton, N.H. Genetics 205, 1335–1351 (2017).
"""
function ϕ(t::Real, r::Real, D_e::Function, σ::Real)
    if iszero(t)
        return zero(eltype(t))
    end
    σ² = σ^2
    denom = 4 * π * t * σ²
    exponent = -r^2 / (4 * t * σ²)
    return (1 / (2 * D_e(t))) * (1 / denom) * exp(exponent)
end

"""
    expected_ibd_blocks_constant_density(r::Real, D::Real, σ::Real, L::Real, G::Real) -> Real

Computes the expected number of identity-by-descent (IBD) blocks of length `L` for a model with constant population density.

```math
\\mathbb{E}[N_L] = \\int_0^\\infty \\mathbb{E}[N_L^t] \\,dt =
\\frac{G}{8\\pi D \\sigma^2}
\\left(\\frac{r}{\\sqrt{L} \\sigma}\\right)^2
K_2\\left(\\frac{\\sqrt{2L} \\, r}{\\sigma}\\right)
```

where:
- `r` is the geographic distance between samples,
- `D` is the effective population density (diploid individuals per unit area),
- `σ` is the root mean square dispersal distance per generation,
- `L` is the minimum length of the IBD block (in Morgans),
- `G` is the total map length of the genome (in Morgans),
- `K₂` is the modified Bessel function of the second kind of order 2.

Reference:
Ringbauer, H., Coop, G., & Barton, N. H. (2017). Genetics, 205(3), 1335–1351.
"""
function expected_ibd_blocks_constant_density(r::Real, D::Real, σ::Real, L::Real, G::Real)
    # Input validation
    if r < 0 || D ≤ 0 || σ ≤ 0 || L ≤ 0 || G ≤ 0
        throw(ArgumentError("All parameters must be positive"))
    end
    term1 = G / (8π * D * σ^2)
    term2 = (r / (sqrt(L) * σ))^2
    term3 = safe_adbesselk(2, sqrt(2L) * r / σ)
    term1 * term2 * term3
end

"""
    expected_ibd_blocks_power_density(r::Real, D::Real, β::Real, σ::Real, L::Real, G::Real) -> Real

Computes the expected number of identity-by-descent (IBD) blocks of length `L` for a model with a power population density in the form of ``D(t) = D_0t^{-β}``

```math
\\mathbb{E}[N_L] = \\int_0^\\infty \\mathbb{E}[N_L^t] \\,dt =
2^{\\frac{-3\\beta}{2}-3}
\\frac{G}{\\pi D \\sigma^2}
\\left(\\frac{r}{\\sqrt{L} \\sigma}\\right)^{(2+\\beta)}
K_{2+\\beta}\\left(\\frac{\\sqrt{2L} \\, r}{\\sigma}\\right)
```

where:
- `r` is the geographic distance between samples,
- `D` is the effective population density (diploid individuals per unit area),
- `β` is the growth rate of the population
- `σ` is the root mean square dispersal distance per generation,
- `L` is the minimum length of the IBD block (in Morgans),
- `G` is the total map length of the genome (in Morgans),
- `K₂` is the modified Bessel function of the second kind of order 2.

Reference:
Ringbauer, H., Coop, G., & Barton, N. H. (2017). Genetics, 205(3), 1335–1351.
"""
function expected_ibd_blocks_power_density(
    r::Real,
    D::Real,
    β::Real,
    σ::Real,
    L::Real,
    G::Real,
)
    # Input validation
    if r < 0 || D ≤ 0 || σ ≤ 0 || L ≤ 0 || G ≤ 0
        throw(ArgumentError("All parameters except β must be positive"))
    end
    term1 = 2^(-3β / 2 - 3)
    term2 = G / (π * D * σ^2)
    term3 = (r / (sqrt(L) * σ))^(2 + β)
    term4 = safe_adbesselk(2 + β, sqrt(2L) * r / σ)
    term1 * term2 * term3 * term4
end

end
