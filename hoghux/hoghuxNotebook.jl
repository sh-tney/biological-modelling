### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ ba8d4db8-a7c0-11eb-3a19-031e97d7051a
md"# 1.1 Intro"

# ╔═╡ 655e3915-11f2-4adf-898e-8fb177b9ea44
md"**Single-Channel Current-Voltage Relationship:** Each channel is either open (conducting) or
closed (non-conducting). Individual channel state transitions are very fast, in the order of 10s of
microseconds. Individual channel conductances are in the order of 10pS.
Voltages are given here in mV referenced to resting membrane potential (RMP), following
Hodgkin and Huxley. RMP in the squid giant axon (internal potential referenced to external
ground) is -72mV."

# ╔═╡ 52d161eb-06a7-404d-9797-d35e3c4575e7
begin
	g_Na = 1.0e-11 	# single sodium channel conductance ~10pS
	E_Na = 115.0 	# sodium reversal potential
	g_K = 1.5e-11 	# single potassium channel conductance ~15pS
	E_K = -12.0 	# potassium reversal potential
	g_L = 2.0e-12 	# single leak channel conductance ~2pS
	E_L = 10.6 		# leak reversal potential
end

# ╔═╡ 95333cb1-c831-4d91-932b-c51262b47fe1


# ╔═╡ Cell order:
# ╟─ba8d4db8-a7c0-11eb-3a19-031e97d7051a
# ╟─655e3915-11f2-4adf-898e-8fb177b9ea44
# ╠═52d161eb-06a7-404d-9797-d35e3c4575e7
# ╠═95333cb1-c831-4d91-932b-c51262b47fe1
