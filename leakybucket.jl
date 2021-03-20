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
md"**The Leaky Bucket Model**"

# ╔═╡ c0f9a8ba-8925-11eb-0def-09b673bc8bee
load("leaky_bucket.png")

# ╔═╡ 7330afdc-8927-11eb-1dc5-63fe6b664e44
md"The leaky bucket model ia a concrete hydrodynamic analogy for how membrane potential depends on currents entering and leaving a cell, and how currents depend on membrane capacitance, channel conductances and voltages. 

Like any model, the leaky bucket model is useful for learning how stuff works, and how to simulate stuff, but it is simplistic and misleading if taken too seriously. 

We will use the leaky bucket model to learn some basic concepts of neuronal function,  and how to use Julia to simulate the underlying physical mechanisms. We can then build on what we learned using this simple model to develop more complex, realistic models of biophysical mechanisms in sensory receptor cells and neurons.

In the leaky bucket model above, the water level represents the membrane potential of a neuron. Current enters at flow rate $u_{IN}$ and leaves at flow rate $u_{OUT}$. The potential rises (or falls) at a rate proportional to the net current new flow:"

# ╔═╡ 0dc99250-8925-11eb-1df0-23c9f390a7a0
md"$\frac{dv}{dt} = \alpha ( u_{IN} - u_{OUT})$"

# ╔═╡ ad49266c-8928-11eb-08e1-5b3f9072c114
md"The out flow rate ($u_{OUT}$) is proportional to the difference between the water level ($v$) and its resting level ($v_{REST}$), and proportional to the cross-sectional area of the leak channel (i.e. how easily the channel conducts current):"

# ╔═╡ 413777de-8929-11eb-191e-0fec3ca38e38
md"$u_{OUT} = \lambda(v - v_{REST})$"

# ╔═╡ Cell order:
# ╠═f195b91c-8909-11eb-0693-552057fb2a76
# ╟─2dc87742-8927-11eb-357f-eb4a35b172d3
# ╟─c0f9a8ba-8925-11eb-0def-09b673bc8bee
# ╟─7330afdc-8927-11eb-1dc5-63fe6b664e44
# ╟─0dc99250-8925-11eb-1df0-23c9f390a7a0
# ╟─ad49266c-8928-11eb-08e1-5b3f9072c114
# ╟─413777de-8929-11eb-191e-0fec3ca38e38
