# # Theory overview
# This package provides an efficient implementation of the inference scheme proposed by Ringbauer et al. ( 2017) to estimate the mean dispersal rate and the effective population density of a population.
# Here, we present an overview of the theory behind the method, but we refer to the original publication for details.

# Mates tend to live close to each other and their offspring. This results in an inverse correlation between geographical distance and genetic relatedness. The goal of this package is to do demographic inference in such spatial genetic patterns.

# ## Diffusion approximation of dispersal

# First, let’s provide some intuition on how we model the dispersal of individuals in a continuous space.
# From an ecological standpoint, we might be interested in the details of the single-generation dispersal, such as modeling the mating process or differences in dispersal across individuals and sexes.
# This is, however, not generally feasible to estimate from genetic data. Instead, we aim to estimate the mean dispersal rate across generations.

# This method approximates the spatial movement of genetic material using a diffusion process.
# Let’s denote as $r(t)$ the distance between two homologous loci at time $t$.
# First, recall that when $t=TMRCA$ (time to the most common recent ancestor), then $r(\text{TMRCA}) = 0$. The distance
# between both loci at the present, $r(0)=r_{\text{obs}}$ (i.e. the distance we can observe), is the sum of the sequence of migration events since the time to the most common recent ancestor.

# If migrations are independent and identically distributed, then the random variable $r_{\text{obs}}$ converges to its own Gaussian distribution with axial variance of $\sigma^2\cdot \text{TMRCA}$ as $\text{TMRCA}$ grows according to the central-limit theorem. Here $sigma^2$ is the average squared axial parent-offspring distance (i.e. the dispersal rate).

# This result is independent of the details of the single-generation dispersal process, while the variance is finite. For example, let’s consider a population where the single-generation dispersal can be modeled as a uniform distribution. Below, we visualize the distribution of distances between homologous loci when $ t = 1$ (single-generation dispersal) and $ t=30$.
# For simplicity, we assume a 1-dimensional space. As we can see, the approximation becomes very accurate even
# after a few generations.

using Distributions, Plots, StatsPlots
displacement = Uniform(-1, 1)
σ2 = var(displacement) # Theoretical variance of displacement
distance(t) = sum(rand(displacement, t)); # Distance after t migrations
n_draws = 10_000
p1 = histogram(
    [distance(1) for _ = 1:n_draws],
    normalize = true,
    label = "t = 1",
    xlabel = "Distance",
    ylabel = "Density",
)
plot!(p1, Normal(0, sqrt(σ2 * 1)), lw = 2, label = "Gaussian Approximation")

p2 = histogram(
    [distance(30) for _ = 1:n_draws],
    normalize = true,
    label = "t = 30",
    xlabel = "Distance",
    ylabel = "Density",
)
plot!(p2, Normal(0, sqrt(σ2 * 30)), lw = 2, label = "Gaussian Approximation")
plot(p1, p2, layout = (1, 2), size = (900, 400), legend = :topright, framestyle = :box)

# This approximation is expected to work well in regimes where the time to the most common recent ancestor is not extremely small
# and when migrations are not correlated and time-homogeneous (i.e. they are independent and identically distributed).

# From the previous section, it is clear that it would be possible to design an inference scheme based on the diffusion approximation
# if we knew the exact time to the most recent common ancestor between homologous loci. For example, the
# probability of two homologous loci being one unit apart given that they coalesced 20 generations ago can
# be calculated from the probability density of the Gaussian approximation for any σ.
pdf(Normal(0, sqrt(σ2 * 20)), 1)

# ## Coalescent theory

# However, the time to the most common recent ancestor is not directly observable. Instead, what we can do is
# marginalize across all possible times $[0, \infty)$. Here is when coalescent theory is used to model the probabilities
# of coalescent events. Roughly speaking, we consider a scenario where two homologous loci *might* coalesce when they
# become very close. The rate with which they coalesce when they are very close is inversely proportional to local effective density $D_e$.
# According to this model, the probability that two homologous loci that are $d$ units apart in the present coalesce at time $t$
# is simply the product of the probability that $r(t) = 0$ and a rate of local coalescence $\frac{1}{D_e(t)}$.

generations = 1:300
D_e = 1 # Constant effective density of 1 (individuals per unit area)
σ = 1 # Constant dispersal rate of 1( distance units per generation)
r = 3 # Pairwise distance between loci
densities = [pdf(Normal(0, sqrt(σ^2 * t)), r) / D_e for t in generations]
plot(
    generations,
    densities,
    ylabel = "Probability of coalescence ϕ(t)",
    xlabel = "Generations",
    label = "(D=1, σ=1, r=3)",
)

# ## IBD blocks
# The last ingredient of this method is the relationship between identity-by-descent blocks
# and time to the most recent common ancestor. First, we point out that generally only recent coalescent
# events are informative of the dispersal rate. This is the reason why this inference scheme relies on
# identity-by-descent blocks, as they can be used to model recent demography.
#
# We say that a segment of DNA is a shared identity-by-descent block if it has been inherited from a common ancestor
# without being broken by recombination. A shared identity-by-descent block of length $l$ that finds its common ancestor at time $t$
# has experienced $2t$ meiosis. Therefore, we expect $l$ to decay with $t$.
#
# More specifically, if we model recombination as a Poisson process and measure genetic distance in Morgans (as it is typically done), then
# the probability that a region of length $ L$ does not recombine follows the exponential distribution with rate
# $\exp(-2Lt)$
#
# Ringbauer et al. (2017) derived that the expected density of identity-by-descent blocks of length $l$ per pair of haploid individuals and
# per unit of block length given that the time to the most recent common ancestor is $t$ is given by
# $E[K_l | t] \approx G 4 t^2 \exp(-2lt)$
# where $G$ is the length of the genome in Morgans. This equation is obtained by accounting for the probability that there is a region of length $l$
# that has not recombined (the $\exp(-2lt)$ term) and it is delimited by two recombination events (the $4 t^2$ term) and summing
# across all possible start sites (the $G$ factor). A slightly more complex expression that accounts for the effect of chromosomal edges
# and diploidy is provided in the article and the software implementation.

# ## Inference scheme
#
# As mentioned above, the inference scheme relies on computing the expected density of identity-by-descent blocks under a given
# demographic model. This involves marginalizing across all possible times of coalescence, therefore computing integrals in the form of
#
# $E[L] = \int_0^\infty E[L | t] \phi(t) dt$
#
# where $f(t)$ is the probability density function of the time to the most recent common ancestor. This expression is hard to track
# down analytically. Ringbauer et al. (2017) derived an analytical solution for a family of demographic models. Alternatively, this
# software implementation provides a numerical approximation of the integral using Gaussian-quadrature rules. We refer to the
# rest of the sections of the documentation for more details.
#
#  Ringbauer et al. (2017) proposed an inference scheme where they assumed the number of observed shared-identity-by-descent blocks
#  whose length fall in a small bin follows a Poisson distribution. The rate of the distribution can be calculated from the expected
#  density of identity-by-descent blocks under the demographic model of interest. If the bin is small enough, it can be simply calculated
#  as $E[L, L+\Delta L ] = E[L] \Delta L $.
#
#  Finally, they propose to approximate the likelihood of the observed data across many bins and pairs of individuals using a composite
#  likelihood (that is, by assuming pairwise observations are independent).
#
#  $L(\theta | Y_i^j) = \Pr(K = Y_i^j | E[L_i]^j \Delta L)$
#
# where $Y_i^j$ is the number of observed shared-identity-by-descent blocks whose length fall in the $i$th bin and are
# shared by the $j$th pair of individuals. The composite likelihood of all data is simply:
#
#  $CL(\theta | \text{Data}) = \prod_{i,j} \Pr(K = Y_i^j | E[L_i]^j \Delta L)$
