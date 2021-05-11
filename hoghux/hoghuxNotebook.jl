### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 30ce572c-660e-4cdb-bcb3-3c40ff7b7d6e
begin
	using PlutoUI
	using Images
	using ImageMagick
	using Plots
	using Unitful
	using UnitfulRecipes
	gr()
end

# ╔═╡ b68f7757-d999-4d42-8133-6a77b2048aba
md"# The Hodgkin-Huxley Model"

# ╔═╡ ba8d4db8-a7c0-11eb-3a19-031e97d7051a
md"## 1.1: Intro"

# ╔═╡ 655e3915-11f2-4adf-898e-8fb177b9ea44
md"**Single-Channel Current-Voltage Relationship:** Each channel is either open (conducting) or
closed (non-conducting). Individual channel state transitions are very fast, in the order of 10s of
microseconds. Individual channel conductances are in the order of 10pS (Siemens are the measure of conductance - inverse of ohms for resistance).
Voltages are given here in mV referenced to resting membrane potential (RMP), following
Hodgkin and Huxley's original data. RMP in the squid giant axon (internal potential referenced to external
ground) is -72mV."

# ╔═╡ 52d161eb-06a7-404d-9797-d35e3c4575e7
begin
	g_Na = 1.0e-11u"S" 			# single sodium channel conductance ~10pS
	E_Na = 115.0u"mV"			# sodium reversal potential
	g_K = 1.5e-11u"S" 			# single potassium channel conductance ~15pS
	E_K = -12.0u"mV"			# potassium reversal potential
	g_L = 2.0e-12u"S"	 		# single leak channel conductance ~2pS
	E_L = 10.6u"mV"				# leak reversal potential
end

# ╔═╡ c729a5ff-b93e-4916-80a8-f9cf717e9ccf
begin
	Vsingle = -50u"mV":1u"mV":150u"mV"	# Voltage array, stepping by 1mV
	
	plot(Vsingle, g_Na*(Vsingle.-E_Na),		# Mapping channel current per Voltage
		label="Na",
		ylabel="mA",
		xticks=(-50:25:150),
		yticks=(-2.5e-9:0.5e-9:3.0e-9),
		legend=:bottomright,
		framestyle=:origin,
		title="Single channel current-voltage plots")
	
	plot!(Vsingle, g_K*(Vsingle.-E_K),
		label="K")
	
	plot!(Vsingle, g_L*(Vsingle.-E_L),
		label="L")
end

# ╔═╡ f093b868-d222-444b-b3f7-dd65efe5642a
md"From the above graph, we can see the differences between channel conductances and their effects on that channels current/voltage relationship. 

When each of these lines cross the x-axis, is the reversal potential (a voltage), at which a current equilibrium is reached where current is neither flowing inward or outward of the channel - the net current is zero. We call the voltage this occurs at  $E$ for each channel.

The equations used to generate these plots all happen to be the components of a simple version of the Hodgkin-Huxley model:

$C\frac{dV}{dt} = -G_{Na}(V - E_{Na}) - G_{K}(V - E_K) - G_L(V - E_L)$

Much like the Leaky Bucket model, this version is simple because it assumes that the conductances are constant - which they are not."

# ╔═╡ ea27dae8-3e0d-42c5-981a-5723d4a1b1e5
md"## 1.2: Channel Open-State Probabilities"

# ╔═╡ ee6bb5e5-d570-48c6-9603-ce224c4e0018
md"Conductance in a whole patch, compartment or neuron is built up from some number of individual voltage-gated channels. Calculating this depends on individual channel conductance ($g_x$) times number of channels ($N$) times the probability that a channel is open at a given voltage ($P_o$):

$G_x(V) = Ng_xP_o(V)$

Single channel conductance and number of channels can be estimated by statistical analysis of conductance in a membrane patch (at a fixed voltage, the total current has a binomial distribution, scaled by channel conductance times voltage). Models are usually parameterized by total conductance density ($mS/cm^2$), assuming all channels are open ($g_x^∗$) times membrane area ($A$):

$\bar{g}_x = Ag^*$

With that in mind, we can substitute this back into our first equation in place of channel count times conductance:

$G_x(V) = \bar{g}_xP_o(V)$

When converting, current densities are in the order of $100mA/cm^2$, or $μA$ for a large cell."

# ╔═╡ 2738ebf4-9433-4f40-8ec5-aaf956098d88
md"### 1.2.1: Equilibrium Open-State Probability"

# ╔═╡ a0e8a0e8-5a67-44e4-990e-924ecf275d97
md"Single-channel open-closed transitions are effectively instantaneous, but they occur probailistically when thermal noise pushes a channel from one state to the other over the energy barrier, and it takes a while (on the timescale of neural computation) to reach equilibrium. The transition rate is faster from the high- to the low-energy state. Equilibrium is reached when the proportion of channels in each state matches the ratio of rate constants for transitions to that state. This typically occurs with a time constant in the order of milliseconds."

# ╔═╡ e5f4b0e7-dbce-496c-b9c1-6bca36011d0b
load("boltzmann.png")

# ╔═╡ 4e4fbb13-6ce0-497a-8ae9-550d3410132f
md"Earlier e assumed that the open state probability is always equal to the equilibrium probability at a given voltage ($P_∞(v)$), because we assumed
that changes on a 1ms timescale are irrelevant. Substututing $P_C = 1 − P_O$
into Boltzmann’s equation for the relative probability of being in one of two states at thermal equilibrium (above) gives a sigmoidal relationship between stimulus energy and open probability at equilibrium:

$P_∞(ΔE) = \frac{1}{1+e^{ΔE/K_BT}}$"

# ╔═╡ 57205af2-edf9-4329-a7ba-ddcf1b8fc540
md"### 1.2.2: Potassium Channel Kinetics"

# ╔═╡ a106af42-ac3d-4241-9161-cb70aa14fd78
md"However, because action potentials are typically about 1ms across, the kinetics of channel opening are essential in modelling action potentials. Using $P_C = 1 − P_O$, we can write a single ODE for channel opening kinetics:

$\frac{dP_o}{dt} = r_o(1 - P_o) - r_cP_o$

Hodgkin and Huxley discovered voltage-dependent sodium and potassium conductances in
squid giant axons at Plymouth Marine Lab in 1939, but their work was interrupted by World War II and was not published until 1952. Advances in electronics and control theory during the war made voltage clamping possible, so that by 1950 they were able to quantify the relationship between trans-membrane voltages and conductances for sodium and potassium ions.

Ion channels were unknown at the time. Hodgkin and Huxley had initially assumed that ions must be transported across the membrane in lipid carriers. However, they noted that the empirically observed kinetics of conductance change could be explained by a Boltzmann-type thermodynamic model, and therefore hypothesised that unspecified “particles” in the membrane have conducting and non-conducting states with different energy levels using:

- m to represent sodium channel open state probability 
- n to represent potassium open state probability 
- α to represent opening rate constant 
- β to represent the closing rate constant

Thus their ODE for the open state probability of potassium conductance “particles” is:

$\frac{dn}{dt} = α_n(1-n) - β_nn$

This predicts an exponential-decay trajectory towards a new equilibrium probability, $P_∞(v)$, following a step change in trans-membrane voltage (like in the Leaky Bucket Model). But this model did not fit the voltage-clamp data. 

By laborious calculation using a mechanical calculator to solve ODES (there
was one programmable electronic computer in Britain at the time, at Manchester University) they showed that membrane conductance is proportional to $n^4$. The current through voltage-gated potassium channels is given by:

$I_K = -\bar{g}_Kn^4(V - E_K)$

where $\bar{g}_K$ is the total conductance if all channels are open (from earlier). From the $n^4$, Hodgkin and Huxley deduced that
“particles” must come in sets of independantly probabalistc channels, and they must all jump into the conducting state in order for current to flow. Decades later it was discovered that voltage-gated channels are tetramers and that indeed all four voltage-sensing elements must be activated for current to flow through the channel."

# ╔═╡ bae9e2f7-c6d2-4674-bd21-be1a12d58401
md"### 1.2.3: Sodium Channel Kinetics"

# ╔═╡ 6bc109d2-e202-4a53-bf7a-2848b2640af6
md"""Voltage-gated sodium channels turn out to be slightly more complicated. Hodgkin and Huxley had to use two types of gating “particle” in order to fit their data. The model requires three particles of one type and one of the other type to be in the conducting state.

$I_{Na} = -\bar{g}_{Na}m^3h(V - E_{Na})$

The second type of particle ($h$) was called “inactivating” because it works backwards, it closes when the membrane is depolarized rather than to open like the rest ($m$).

We now know that voltage-gated sodium channels, like voltage-gated potassium channels, have four subunits with independently triggered voltage sensors, and a single, independent, inactivation gate that blocks the channel.

In that case you might expect that the exponent of $m$ should be 4, not 3, in the Hodgkin-Huxley model. The short answer is that you can get realistic-looking simulated action potentials using an exponent of 4, but then it wouldn’t be the Hodgkin-Huxley model, would it? 

The longer answer is that the reality is more complicated, and increasing the exponent of m from 3 to 4 in the Hodgkin-Huxley model doesn’t make it usefully more accurate or realistic.

If we have a question whose answer might depend subtley on the precise timing and shape of action potentials, we will need to construct a more realistic model. The Hodgkin-Huxley model is a good foundation for learning how to do that. Biophysical models of neurons tend to be called “Hodgkin-Huxley models” even if they are not actually the Hodgkin-Huxley model."""

# ╔═╡ 62fb9bbf-1273-4cd3-9950-b83c50b4625b
md"### 1.2.4: Leak Channels"

# ╔═╡ 81bd104e-ab89-4dcc-b902-64ac0bd8f8b7
md"Hodgkin and Huxley also reported a leak current, due to a conductance that does not depend on voltage:

$I_L = -\bar{g}_L(V - E_L)$"

# ╔═╡ 857627b0-f86c-4756-991d-b89c71edbd7c
md"## 1.3: The Hodgkin Huxley Model"

# ╔═╡ 4306dbf4-ebe9-43d9-bf88-6a4b51c50d39
md"$C\frac{dV}{dt} = -\bar{g}_{Na}m^3h(V-E_{Na}) - \bar{g}_Kn^4(V-E_K) - \bar{g}_L(V-E_L)$

$\frac{dm}{dt} = α_m(1-m) - β_mm$

$\frac{dn}{dt} = α_n(1-n) - β_nn$

$\frac{dh}{dt} = α_h(1-h) - β_hh$

This is a nonlinear dynamic model in which membrane potential affects current flow, which in turn, affects membrane potential. There is an equilibrium when the net current (the sum of all terms on the right) is zero. 

The voltage at equilibrium is called the resting membrane potential. Because
the sodium equilibrium potential is above resting potential, voltage-gated sodium channels cause positive feedback; depolarization increases the depolarizing current.

The open state probability for each channel has an equilibrium when the right hand side adds to zero, for example:

$m_∞(V) = \frac{α_m}{α_m + β_m}$

Each of these approach the equilibrium at time constant $τ$:

$τ_m(V) = \frac{1}{α_m + β_m}$

These describe exponential decay towards m∞(v) with time constant τm(v). The ODE below allows us to model dynamic changes in conductance while the membrane potential is changing (including feedback between membrane potential and conductance):

$τ_m(V)\frac{dm}{dt} = m_∞(V) - m$"

# ╔═╡ ec110792-cfa4-4b99-b2ae-204ad77ca37b
md"### 1.3.1: Activation Curves"

# ╔═╡ 79f2588d-c23c-490c-b13c-2921ad0624a2
md"The activation/inactivation curve for a channel is a plot of the equilibrium open-state probability,
$P_∞(V)$ for activation/inactivation gates."

# ╔═╡ 71f8325a-8ba0-4656-8f61-919ea1699400
begin
	# α: Opening rate constant
	# β: Closing rate constant
	# τ: Time constant to approach equilibrium
	
	# All terms of v are divided by "mV" to cancel out the units so Unitful
	# works, and because these are probabilities they don't matter anyway.
	
	# n: Potassium channel
	α_n(v) = v/u"mV" == 10.0 ? .1 : 0.01*(10.0-v/u"mV")/(exp((10.0-v/u"mV")/10.0)-1.0)
	β_n(v) = 0.125*exp(-v/80.0u"mV")
	τ_n(v) = 1.0u"ms"/(α_n(v) + β_n(v))
	n_inf(v) = α_n(v)*τ_n(v)/u"ms"
	
	
	# m: Sodium activation channel
	α_m(v) = v/u"mV" == 25.0 ? 1.0 : 0.1*(25.0-v/u"mV")/(exp((25.0-v/u"mV")/10)-1.0)
	β_m(v) = 4.0*exp(-v/18.0u"mV")
	τ_m(v) = 1.0u"ms"/(α_m(v) + β_m(v))
	m_inf(v) = α_m(v)*τ_m(v)/u"ms"
	
	
	# h: Sodium inactivation channel
	α_h(v) = 0.07*exp(-v/20.0u"mV")
	β_h(v) = 1.0/(exp((30.0-v/u"mV")/10)+1.0)
	τ_h(v) = 1.0u"ms"/(α_h(v) + β_h(v))
	h_inf(v) = α_h(v)*τ_h(v)/u"ms"
end

# ╔═╡ f9e94324-4a7e-4d28-be75-2c022379e6bd
begin
	Vp = -50u"mV":1u"mV":150u"mV"	# Voltage array, stepping by 1mV
	
	plot(Vp, n_inf.(Vp),			# Mapping channel open probability per Voltage
		label="K Activation",
		ylabel="Probability",
		legend=:right,
		title="Channel sub-unit open-state Probabilities")
	
	plot!(Vp, m_inf.(Vp),
		label="Na Activation")
	
	plot!(Vp, h_inf.(Vp),
		label="Na Inactivation")
end

# ╔═╡ 08379b09-db4a-479b-bf0b-d5e5636ac5bc
begin
	Vτ = -50u"mV":1u"mV":150u"mV"	# Voltage array, stepping by 1mV
	
	plot(Vτ, τ_n.(Vτ),				# Mapping equilibrium time constant 
		label="K Activation",
		ylabel="τ",
		legend=:right,
		title="Channel sub-unit to-Equilibrium Time Constants ")
	
	plot!(Vτ, τ_m.(Vτ),
		label="Na Activation")
	
	plot!(Vτ, τ_h.(Vτ),
		label="Na Inactivation")
end

# ╔═╡ 66219bb1-eae7-45cb-ad4f-f48b72423366
md"### 1.3.2: Current-Voltage Plots"

# ╔═╡ e3b3f228-dcbb-43d5-9077-a95885d4ac0b
begin
	# Membrane conductance constants, from the original Hodgkin-Huxley paper
	gstar_Na = 120.0u"mS/cm^2"
	gstar_K = 36.0u"mS/cm^2"
	gstar_L = 0.3u"mS/cm^2"

	d = 50.0e-3u"cm"	# Cell diameter 500um (This is a squid giant neuron)
	A = π*d^2 			# Membrane area
	
	# Cell conductance in mS
	gbar_Na = A*gstar_Na
	gbar_K = A*gstar_K
	gbar_L = A*gstar_L
end

# ╔═╡ 50877c70-e7b0-4e61-b03f-9f67c989bed4
begin
	Vg = -50u"mV":1u"mV":150u"mV"	# Voltage array, stepping by 1mV
	
	plot(Vg, gbar_L*(Vg.-E_L),		
		label="Leak",
		ylabel="μA",
		legend=:topleft,
		title="Channel sub-unit current-voltage plots ")
	
	plot!(Vg, gbar_K*(Vg.-E_K).*(n_inf.(Vg)),
		label="K")
	
	plot!(Vg, gbar_Na*(Vg.-E_Na).*(m_inf.(Vg)),
		label="Na")
end

# ╔═╡ c86f9ae6-b593-4512-8fc7-a98d910d338b
md"## 1.4: Numerical Solution of Hodgkin-Huxley"

# ╔═╡ cd9d9afb-6586-4dd9-8ee5-59e422a460c1
md"The state variables of the Hodgkin-Huxley model are v, m, h and n. The state vector is $x = [v, m, h, n]$ (Any order, as long as we keep track of which state variable
correponds to which parameter in the model). The state update function or ODE function returns the derivative (rate of change) of state as a function of state and external input(s). (Always possible for any set of ODEs consistent with causality).

Below is a state-step function, along with a slider to change in the input current for the final output graph."

# ╔═╡ f558b5a9-c5de-40a3-ad20-28af343c739f
# state vector x = [v, m, h, n]
function hh_state_step(x, c, Δt, u_in)
	v = x[1]
	m = x[2]
	h = x[3]
	n = x[4]
	
	Δx = x
	Δx[1] = (-gbar_Na*m^3*h*(v-E_Na) -gbar_K*n^4*(v-E_K) -gbar_L*(v-E_L) +u_in)/c*Δt
	Δx[2] = (m_inf(v) - m)/τ_m(v)*Δt
	Δx[3] = (h_inf(v) - h)/τ_h(v)*Δt
	Δx[4] = (n_inf(v) - n)/τ_n(v)*Δt
	return Δx
end

# ╔═╡ 2e556f89-7ccd-48e5-818a-5c7b789a4735
@bind u Slider(0.0:0.05:3.0, default=0.1)

# ╔═╡ a544b61e-d394-4b9f-a375-2e49d04c304c
begin
	# specific capacitance 1uF/cm^2
	Cs = 1.0u"μF/cm^2"
	C = A*Cs
	
	Δt = 0.01u"ms"
	T = 50.0u"ms"
	t = 0u"ms":Δt:T
	u_in = u * u"mS*mV" # Input current (μA)
	
	# Establish state vectors for all of the changing variables
	n = length(t)
	V = zeros(n)u"mV"
	M = zeros(n)
	H = zeros(n)
	N = zeros(n)
	
	# state vector x = [v, m, h, n]
	for i in 2:n
		X = hh_state_step([V[i-1], M[i-1], H[i-1], N[i-1]], C, Δt, u_in)
		V[i] = V[i-1] + X[1]
		M[i] = M[i-1] + X[2]
		H[i] = H[i-1] + X[3]
		N[i] = N[i-1] + X[4]
	end
	
	plot(t, V, 
		label=false,
		title="Neuron Voltage, Input: $u μA")
end

# ╔═╡ Cell order:
# ╟─b68f7757-d999-4d42-8133-6a77b2048aba
# ╠═30ce572c-660e-4cdb-bcb3-3c40ff7b7d6e
# ╟─ba8d4db8-a7c0-11eb-3a19-031e97d7051a
# ╟─655e3915-11f2-4adf-898e-8fb177b9ea44
# ╠═52d161eb-06a7-404d-9797-d35e3c4575e7
# ╠═c729a5ff-b93e-4916-80a8-f9cf717e9ccf
# ╟─f093b868-d222-444b-b3f7-dd65efe5642a
# ╟─ea27dae8-3e0d-42c5-981a-5723d4a1b1e5
# ╟─ee6bb5e5-d570-48c6-9603-ce224c4e0018
# ╟─2738ebf4-9433-4f40-8ec5-aaf956098d88
# ╟─a0e8a0e8-5a67-44e4-990e-924ecf275d97
# ╟─e5f4b0e7-dbce-496c-b9c1-6bca36011d0b
# ╟─4e4fbb13-6ce0-497a-8ae9-550d3410132f
# ╟─57205af2-edf9-4329-a7ba-ddcf1b8fc540
# ╟─a106af42-ac3d-4241-9161-cb70aa14fd78
# ╟─bae9e2f7-c6d2-4674-bd21-be1a12d58401
# ╟─6bc109d2-e202-4a53-bf7a-2848b2640af6
# ╟─62fb9bbf-1273-4cd3-9950-b83c50b4625b
# ╟─81bd104e-ab89-4dcc-b902-64ac0bd8f8b7
# ╟─857627b0-f86c-4756-991d-b89c71edbd7c
# ╟─4306dbf4-ebe9-43d9-bf88-6a4b51c50d39
# ╟─ec110792-cfa4-4b99-b2ae-204ad77ca37b
# ╟─79f2588d-c23c-490c-b13c-2921ad0624a2
# ╠═71f8325a-8ba0-4656-8f61-919ea1699400
# ╠═f9e94324-4a7e-4d28-be75-2c022379e6bd
# ╠═08379b09-db4a-479b-bf0b-d5e5636ac5bc
# ╟─66219bb1-eae7-45cb-ad4f-f48b72423366
# ╠═e3b3f228-dcbb-43d5-9077-a95885d4ac0b
# ╠═50877c70-e7b0-4e61-b03f-9f67c989bed4
# ╟─c86f9ae6-b593-4512-8fc7-a98d910d338b
# ╟─cd9d9afb-6586-4dd9-8ee5-59e422a460c1
# ╠═f558b5a9-c5de-40a3-ad20-28af343c739f
# ╟─2e556f89-7ccd-48e5-818a-5c7b789a4735
# ╠═a544b61e-d394-4b9f-a375-2e49d04c304c
