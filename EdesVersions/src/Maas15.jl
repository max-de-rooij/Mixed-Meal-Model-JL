using ComponentArrays, SysBioModels

# EDeS Parameters
parameters = ModelParameters(
    estimated = ComponentVector(k1 = 1.35e-2, k5 = 3.80e-3, k6 = 5.82e-1),
    fixed = ComponentVector(k2 = 0.28, k3 = 6.07e-3, k4 = 2.35e-4, k7 = 1.15, k8 = 7.27, k9 = 3.83e-2, k10 = 2.84e-1, sigma = 1.4, Km = 13.2, Gb = 5.0, Ib = 8.0, Ghist = 5.0),
    constants = ComponentVector(f = 0.00551, vg = 17/70, gbliv = 0.043, beta = 1.0, taui = 31.0, taud = 3.0, vi = 13/70, Gthpl = 9., t_integralwindow = 30., c1 = 0.1)
)

# EDeS Equations

# Gut glucose
function dGgut(u, p, t, input)
    @. p.sigma * p.k1^p.sigma * t^(p.sigma - 1) * exp(-(p.k1 * t)^p.sigma) * input[1] - p.k2 * u[1]
end

# Plasma glucose
function dGpl(u, p, input)
    @. (p.gbliv - p.k3 * (u[2] - p.Gb) - p.k4 * p.beta * u[5]) + (p.k2 * (p.f / (p.vg * input[2])) * u[1]) - (p.gbliv * ((p.Km + p.Gb)/p.Gb) * (u[2] / (p.Km + u[2]))) - (p.k5 * p.beta * u[5] * (u[2] / (p.Km + u[2]))) - ((p.c1 / (p.vg * input[2])) * (u[2] - p.Gthpl) * (u[2] > p.Gthpl))
end

# glucose integrator
function dGint(u, p)
    @. u[2]-p.Ghist
end

# Plasma insulin
function dIpl(u, p, duG)
    @. (p.beta^-1) * (p.k6 * (u[2] - p.Gb)  + (p.k7 / p.taui) * (u[3]+p.Gb) + (p.k8 * p.taud) * duG) - (p.k7 * p.Gb / (p.beta * p.taui * p.Ib)) * u[4] - p.k9 * (u[4] - p.Ib)
end

# Remote insulin
function dIrem(u, p)
    @. p.k9*(u[4] - p.Ib) - p.k10*u[5]
end

# Constructor equation
function edes_signature(D::Dict{Symbol, Function}, inputs::T) where T<:AbstractVector 
    sys! = function(du, u, p, t)
        du[1] = D[:GutGlucose](u,p,t,inputs)
        du[2] = D[:PlasmaGlucose](u,p,inputs)
        du[3] = D[:GlucoseIntegrator](u, p)
        du[4] = D[:PlasmaInsulin](u, p, du[2])
        du[5] = D[:RemoteInsulin](u, p)
    end 
    sys!
end

# Insulin Integrator Callback
condition = function (u, t, integrator)
    t > integrator.p.t_integralwindow
end
integratorfun!(integrator) = integrator.p.Ghist = integrator.u[2];

# Define SimulationSpecs
simulation_specs = SimulationSpecs(
    initialvalues = (p, t0) -> [0., p.Gb, 0., p.Ib, 0.],
    timespan = (0., 120.)
)

# Define System
system_specs = SystemSpecs(
    Dict(
        :GutGlucose         =>  dGgut,
        :PlasmaGlucose      =>  dGpl,
        :GlucoseIntegrator  =>  dGint,
        :PlasmaInsulin      =>  dIpl,
        :RemoteInsulin      =>  dIrem),
    edes_signature
)

# Define Model
EDES_Maas15 = SysBioModel(
    system_specs,
    simulation_specs,
    parameters
)
