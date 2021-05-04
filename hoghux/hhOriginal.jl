using Unitful
using Plots
using UnitfulRecipes

# membrane conductance consts in mS/cm^2
gstar_Na = 120.0#u"ms/cm^2"
gstar_K = 36.0#u"ms/cm^2"
gstar_L = 0.3#u"ms/cm^2"

# cell diameter 500um (nb insanely big! This is squid giant neuron)
d = 50.0e-3#u"cm" # in cm
A = π*d^2# membrane area
# cell conductance in mS
gbar_Na = A*gstar_Na
gbar_K = A*gstar_K
gbar_L = A*gstar_L

g_Na = 1.0e-11 # single sodium channel conductance ~10pS
E_Na = 115.0 # sodium reversal potential
g_K = 1.5e-11 # single potassium channel conductance ~15pS
E_K = -12.0 # potassium reversal potential
g_L = 2.0e-12 # single leak channel conductance ~2pS
E_L = 10.6 # leak reversal potential

struct Neuron

    #C::Unitful.Quantity{Float64, typeof(u"cm")}
    C::Float64
    V::Float64 # Voltage
    n::Float64
    m::Float64
    h::Float64

end

function Neuron()
    return Neuron(1.0, 0, 0, 0, 0)
end

function Neuron(C)
    return Neuron(C, 0., n_inf(0.), m_inf(0.), h_inf(0.))
end

#formulas from handout
α_n(v) = v == 10.0 ? .1 : 0.01 * (10.0 - v) / (exp((10.0 - v) / 10.0) - 1.0)
β_n(v) = 0.125*exp(-v/80.)
τ_n(v) = 1.0/(α_n(v) + β_n(v))
n_inf(v) = α_n(v)*τ_n(v)
α_m(v) = v == 25.0 ? 1.0 : 0.1*(25.0 - v)/(exp((25.0-v)/10)-1.0)
β_m(v) = 4.0*exp(-v/18.)
τ_m(v) = 1.0/(α_m(v) + β_m(v))
m_inf(v) = α_m(v)*τ_m(v)
α_h(v) = 0.07*exp(-v/20.)
β_h(v) = 1.0/(exp((30.0-v)/10)+1.0)
τ_h(v) = 1.0/(α_h(v) + β_h(v))
h_inf(v) = α_h(v)*τ_h(v)

function HH_update(neuron::Neuron, u::Float64, dt::Float64)

    V = neuron.V
    n = neuron.n
    m = neuron.m
    h = neuron.h
    C = neuron.C

    ΔV = (-gbar_Na*m^3*h*(V - E_Na) - gbar_K*n^4*(V - E_K) - gbar_L*(V - E_L) + u) / C * dt
    Δn = (n_inf(V) - n)/τ_n(V)*dt
    Δm = (m_inf(V) - m)/τ_m(V)*dt
    Δh = (h_inf(V) - h)/τ_h(V)*dt

    return(Neuron(C, V + ΔV, n + Δn, m + Δm, h + Δh))

end

print("hello")

const Δt = 0.01#u"ms"
const T = 50.0#u"ms"
const t = 0:Δt:T #0u"s":Δt:T

n = Neuron(A*1) # 1 F per capacitance

N = length(t)
V = zeros(N)

print("looping")
for i in 1:N
    global n = HH_update(n, 0.1, Δt)
    V[i] = n.V
    println(n.V)
end
print("plotting")
p = plot(t, V)
ylims!(-25, 150)
display(p)
print("bye")

