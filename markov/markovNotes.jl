### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ bc190ffc-b2c9-11eb-2ad1-798abf5245ca
md"When we try to move and simulate the actual expected numbers of subunits opening and closing, we must take a random sample instead of just eyeballing the probability - we take the trajectory of the channel states, not the trajectry of the average across channel states

So instead for each state step of drawing an independant probability of channel state, we will be propogating forward the actual number of open and closed channels, and then evaluating the probability of each opening or closing, and then pushing those states forward to the next step.

$N_{K Closing} = binomial(N_K, β_nΔt)$

$N_{K Opening} = binomial(1-N_K, α_nΔt)$

where $N_K$ is currently open units, and $N-N_K$ are closed units

$\frac{dN}{dt} = α_N(1-N) + β_NN$

$ΔN = α_NΔt(1-N) + β_NΔtN$

When we multiply the probabilities (α & β, which are naturally 1/time-to-switchh) by the time units (Δt), we get the occuring chance that a channel has switched in that time, and we can multiply that by N


Nyquist - Sampling Theory
Euler-Maruyama - Distribution Sample Scaling"

# ╔═╡ Cell order:
# ╠═bc190ffc-b2c9-11eb-2ad1-798abf5245ca
