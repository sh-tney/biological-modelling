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
	using UnitfulRecipes
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

# ╔═╡ 84be2c90-8934-11eb-1fe8-955ab7e8db97
md"Dimensional analysis is a way to check the correctness of a model. The dimensions of all terms in any
model should be the same, and all terms should correspond to real physical quantities."

# ╔═╡ 0f158f3c-8935-11eb-194c-e567a1f60c64
md"$C\frac{dv}{dt} + \lambda(v - v_{REST}) = u_{IN}$"

# ╔═╡ 1da56b44-8935-11eb-1a5f-99ad4536db91
md"""In the leaky bucket model we see that the term on the right hand side, the input current ($u_{IN}$), has dimensions of volume, or length cubed, over time. Therefore the two terms on the left hand side must also have dimensions of length cubed over time. 

In the first term, $C$ has dimensions of length squared, because $\frac{dv}{dt}$ has the dimensions of length over time.

In the second term, $\lambda$ has dimensions of length squared over time ($cm^2.s^{-1}$), because has $v$ the dimension of length ($cm$).

If we choose units of $cm$ (length) and $seconds$ (time) then the $u_{IN}$ is specified in $cm^3$ per second, $\lambda$ is in $cm^2$ per second ($cm^2.s^{-1}$), and $C$ is in $cm^2$. If you look again at the figure, you'll notice that $C$ is the area of the base of the bucket (it is indeed an area) and it makes sense that the leak current is
proportional to the cross-sectional area of the leak channel times the pressure at the mouth of the channel. 

So even without being hydrodynamic experts, we know the dimensions of the model parameters, and have a pretty good idea how to determine plausible values of those  parameters. A dimensionally correct model is not necessarily true, but a model that is not dimensionally correct is not even wrong. It cannot be a realistic description of a physical process.

A model may be dimensionally correct but parameterized inconsistently, e.g. $C$ might be measured in $cm^2$, but $\lambda$ in $\mu m^2.s^{-1}$. Such errors can produce bad simulation results that "look OK", especially if we are not sure what the result should look like. The Julia package *Unitful* can be used to ensure that equations are dimensionally correct and to convert all parameters into consistent units.

The hydrodynamic parameters of the leaky bucket model can be translated to biophysical parameters of a neuron. Water height translates to membrane potential and current (water flow) translates to current (electrical charge flow). The parameter $C$ specifies the rate of water level rise for a given input current, translates to the *membrane capacitance*, which specifies how fast the membrane potential rises for a given input current; it is proportional to the area of the cell membrane. 

The parameter $\lambda$ specifies the current flow in proportion to the difference between the water level and the level of the channel, translates to channel conductance, which specifies the current through the channel in proportion to the membrane potential, referenced to the equiliibrium potential or reversal potential of the channel.

The leaky bucket model helps to develop intuition about neuronal biophysics (up to a point), and is useful for learning how to simulate neuronal biophysics in a more familiar and concrete context."""

# ╔═╡ 0c95b5fe-89ec-11eb-36a4-d7fa11009528
md"## Numerical Solution"

# ╔═╡ 163a4034-89ec-11eb-0da0-738f98d63eed
md"We want to compute the water level over time as we change the input current, and explore how this depends on the model parameters. When we step up from a bucket to a neuron, the question will be how does membrane potential change when the neuron receives a synaptic input current, and how does this depend on the membrane capacitance and channel conductances of the neuron?

To solve the Leaky Bucket Model, we have to first rearrange to have $\frac{dv}{dt}$ on
 the left hand side:"

# ╔═╡ 79682d06-89ec-11eb-28c3-c75547cc182a
md"$\frac{dv}{dt} = (u_{IN} - \lambda(v - v_{REST}))/C$"

# ╔═╡ 8127df7c-89ed-11eb-0508-7b99fda17907
md"This is a simple special case of a *state space model*, in which only first derivatives appear, and each first derivative appears alone on the left hand side of its own equation. It is generally possible to rewrite any system of differential equations as a state space model. Numerical differential equation solvers generally require the equations to be given in this form.

We construct a *finite difference approximation* to the differential equation by replacing the derivative $\frac{dv}{dt}$ with the ratio $\frac{\Delta v}{\Delta t}$, where $\Delta v$ represents a small change in $v$ over a small time interval $\Delta t$:"

# ╔═╡ 5212cfcc-89ee-11eb-16cb-07877e03777b
md"$\frac{\Delta v}{\Delta t} = (u_{IN} - \lambda(v - v_{REST}))/C$"

# ╔═╡ 877b9172-89ee-11eb-1825-e902fe36a505
md"This equation tells us that the change in water level $\Delta v$ over a small time interval $\Delta t$ is:"

# ╔═╡ a043681c-89ee-11eb-2f29-5b679d2ee2c1
md"$\Delta v = (u_{IN} - \lambda(v - v_{REST}))\Delta t/C$"

# ╔═╡ 931e7dec-89ef-11eb-19aa-c3c2bab8eb2e
md"Therefore if we know the water level at some time $t$ we can compute the water level at time $t + \Delta t$:"

# ╔═╡ cc0c0d7c-89ef-11eb-255f-5fec3337bc94
md"$v(t + \Delta t) = v(t) + \Delta v = (u_{IN} - \lambda(v - v_{REST}))\Delta t/C$"

# ╔═╡ 55f3c78c-89f0-11eb-20c1-2763d779cf2e
md"We'll write a Julia function to compute $v(t + \Delta t)$ given $v(t)$. By doing this we can set up a general template for simulating with state space models that stays simple no matter how complicated the differential equation model gets, which is basically why numerical differential equation solvers ask for state space models - once you've done the hard work of writing the model in this form, the rest is easy:"

# ╔═╡ befe34d8-89f0-11eb-27ba-499d66ef9afa
function bucket_state_step(v, Δt, C, λ, v_rest, u)
	Δv = (u - λ*(v - v_rest))*Δt/C
	return Δv
end

# ╔═╡ 7056b634-89f2-11eb-3df1-e5fd49402a93
md"We might guess what will happen to the water level (at least qualitatively) if we
fill the bucket and let the water leak out without adding any more water. 

Let's say that the bucket is cylindrical with radius $15cm$ and it is filled to a height of $11cm$. The leak channel is $1cm$ above the bottom. Its conductance depends on its mouth area and its length, i.e. how hard it is for water to flow along the channel. 

We could do some back of the envelope calculations to estimate a plausible value for $\lambda$, but let's just take a stab at it. I guess that the leak is going to be around $50cm^3$ per second when the bucket is full $(v - v_{REST} = 10)$, which gives $\lambda = 5cm^2.s^{-1}$:"

# ╔═╡ 5d0c4582-89f3-11eb-1e2b-21d39cb751fd
begin
	C = (π*15.0^2)u"cm^2" 		# Cross-sectional area of the bucket
	λ = 5.0u"cm^2/s"			# Coefficient to represent leak flow rate
	Δt = 0.05u"s"				# Simulation running in 5ms steps
	T = 600.0u"s"				# Total duration of simulation - 10 minutes
	t = 0u"s":Δt:T;				# Time iterator
	v = zeros(length(t))u"cm"	# Vector container for plotting computed values
	v[1] = 11.0u"cm"			# Intial water height
	v_rest = 1.0u"cm"			# Leak height
	u_in = 0.0u"cm^3/s"			# Input Current
	v
end

# ╔═╡ 805ee830-89f5-11eb-0051-ed5933e91fa4
begin
	for i in 2:length(t)
		v[i] = v[i-1] + bucket_state_step(v[i-1], Δt, C, λ, v_rest, u_in)
	end
	plot(t, v, xlabel="Time t", ylabel="Water Level v", legend=false, title="Leaky Bucket")
end

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
# ╟─84be2c90-8934-11eb-1fe8-955ab7e8db97
# ╟─0f158f3c-8935-11eb-194c-e567a1f60c64
# ╟─1da56b44-8935-11eb-1a5f-99ad4536db91
# ╟─0c95b5fe-89ec-11eb-36a4-d7fa11009528
# ╟─163a4034-89ec-11eb-0da0-738f98d63eed
# ╟─79682d06-89ec-11eb-28c3-c75547cc182a
# ╟─8127df7c-89ed-11eb-0508-7b99fda17907
# ╟─5212cfcc-89ee-11eb-16cb-07877e03777b
# ╟─877b9172-89ee-11eb-1825-e902fe36a505
# ╟─a043681c-89ee-11eb-2f29-5b679d2ee2c1
# ╟─931e7dec-89ef-11eb-19aa-c3c2bab8eb2e
# ╟─cc0c0d7c-89ef-11eb-255f-5fec3337bc94
# ╟─55f3c78c-89f0-11eb-20c1-2763d779cf2e
# ╠═befe34d8-89f0-11eb-27ba-499d66ef9afa
# ╟─7056b634-89f2-11eb-3df1-e5fd49402a93
# ╠═5d0c4582-89f3-11eb-1e2b-21d39cb751fd
# ╠═805ee830-89f5-11eb-0051-ed5933e91fa4
