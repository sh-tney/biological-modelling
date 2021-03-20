### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ f195b91c-8909-11eb-0693-552057fb2a76
begin
	using PlutoUI
	using Images
	using ImageMagick
	using Plots
	using Unitful
	gr()
end

# ╔═╡ 2dc87742-8927-11eb-357f-eb4a35b172d3
md"# The Leaky Bucket Model"

# ╔═╡ c0f9a8ba-8925-11eb-0def-09b673bc8bee
load("leaky_bucket.png")

# ╔═╡ 7330afdc-8927-11eb-1dc5-63fe6b664e44
md"The leaky bucket model is a concrete hydrodynamic analogy for how membrane potential depends on currents entering and leaving a cell, and how currents depend on membrane capacitance, channel conductances and voltages. 

Like any model, the leaky bucket model is useful for learning how stuff works, and how to simulate stuff, but it is simplistic and misleading if taken too seriously. 

We use the leaky bucket model to learn some basic concepts of neuronal function,  and how to use Julia to simulate the underlying physical mechanisms. We can then build on  this simple model to develop more complex, realistic models of biophysical mechanisms in sensory receptor cells and neurons.

In the leaky bucket model above, the water level represents the membrane potential of a neuron. Current enters at flow rate $u_{IN}$ and leaves at flow rate $u_{OUT}$. The potential rises (or falls) at a rate proportional to the net current flow:"

# ╔═╡ 0dc99250-8925-11eb-1df0-23c9f390a7a0
md"$\frac{dv}{dt} = \alpha ( u_{IN} - u_{OUT})$"

# ╔═╡ ad49266c-8928-11eb-08e1-5b3f9072c114
md"The out flow rate ($u_{OUT}$) is proportional to the difference between the water level ($v$) and its resting level ($v_{REST}$), and proportional to the cross-sectional area of the leak channel (i.e. how easily the channel conducts current):"

# ╔═╡ 413777de-8929-11eb-191e-0fec3ca38e38
md"$u_{OUT} = \lambda(v - v_{REST})$"

# ╔═╡ c98f56c8-892f-11eb-1011-39ae290db9d1
md"Combining these equations and re-arranging to put the dependent quantities on the left and the independent quantities on the right gives (where $C = 1/\alpha$):"

# ╔═╡ dacbbb66-892f-11eb-0912-0fca799b65b6
md"$C\frac{dv}{dt} + \lambda(v - v_{REST}) = u_{IN}$"

# ╔═╡ c19d6f62-8930-11eb-0687-b99ac9ba3831
md"""The above equation is the *"Leaky Integrator Model"*, a standard computational neural model in which $v$ is a given membrane potential, and $v_{REST}$ (or $v_0$) is the resting membrane potential. Membrane potential is measured as **Voltage** (a potential difference) relative to a *ground potential*, which is the potential at some arbitrary point.

Modellers often set this ground potential to zero, corresponding to measuring
voltages relative to the internal potential of a cell when it is at rest, which is the simplest option.

Neurophysiologists on the other hand typically measure membrane potential relative to "ground" potential measured by an electrode "at infnity" - far enough away from the neuron being recorded not to be affected by the activity of the neuron. This is the simplest choice *for an experimenter*. However, it is a potential source of error
and confusion in modelling.

We call the above equation the "*Leaky Bucket Model*" rather than a leaky integrator model, to remind us that although real neurons are indeed leaky integrators that transform currents into voltages, this is a rather clumsy approximation to a real neuron."""

# ╔═╡ 43857056-8934-11eb-30ca-8119300e9390
md"## Dimensional Analysis"

# ╔═╡ Cell order:
# ╠═f195b91c-8909-11eb-0693-552057fb2a76
# ╟─2dc87742-8927-11eb-357f-eb4a35b172d3
# ╟─c0f9a8ba-8925-11eb-0def-09b673bc8bee
# ╟─7330afdc-8927-11eb-1dc5-63fe6b664e44
# ╟─0dc99250-8925-11eb-1df0-23c9f390a7a0
# ╟─ad49266c-8928-11eb-08e1-5b3f9072c114
# ╟─413777de-8929-11eb-191e-0fec3ca38e38
# ╟─c98f56c8-892f-11eb-1011-39ae290db9d1
# ╟─dacbbb66-892f-11eb-0912-0fca799b65b6
# ╟─c19d6f62-8930-11eb-0687-b99ac9ba3831
# ╟─43857056-8934-11eb-30ca-8119300e9390
