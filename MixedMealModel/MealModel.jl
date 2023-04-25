# Building the Meal Model from EDES
using EdesVersions, SysBioModels, ComponentArrays, DifferentialEquations, StaticArrays

edes = EDES(:Maas15)

# Meal Model Specific Parameters
parameters_mealmodel = ModelParameters(
    estimated = ComponentVector(k11 = 0.00045, k12 = 0.0713, tauLPL = 208.88, k14 = 0.0163, k16 = 0.0119),
    fixed = ComponentVector(fspill = 30., ATLmax = 0.215, KATL = 0.0385, k13 = 0.0088, k15 = 1e-5, TGb = 1.3, NEFAb = 0.33),
    constants = ComponentVector(fTG = 0.00113, vTG = 0.06)
)

# Combined parameters
parameters = parameters_mealmodel + edes.parameters

# Meal Model Equations
function dNEFAInsulinDelays(u, p)
    @. [
        (3/p.tauLPL).*(u[4] - u[6]),
        (3/p.tauLPL).*(u[6] - u[7]),
        (3/p.tauLPL).*(u[7] - u[8])
    ]
end

function dTGgut(u, p, t, input)
    @. [
        (p.sigma * p.k13^p.sigma * t^(p.sigma - 1) * exp(-(p.k13 * t)^p.sigma) * input[3]) - p.k14 * u[9],
        p.k14 * u[9] - p.k14 * u[10],
        p.k14 * u[10] - p.k14 * u[11]
    ]
end

function dTGpl(u, p, input)
    @. (p.k16 - p.k15*(u[8]-p.Ib)) + (p.k14 * (p.fTG/(p.vTG*input[2]))*u[11]) - (p.k11*u[12]*u[8])
end

function dNEFApl(u, p)
    @. 3*((1/100)*(p.fspill*(p.Ib/u[6])))*p.k11*u[12]*u[8] + p.ATLmax/(1+p.KATL*(u[6])^2) - p.k12*u[13]
end

# Constructor equation
function mealmodel_signature(D::Dict{Symbol, Function}, inputs::T) where T<:AbstractVector
    sys! = function(du, u, p, t)
        # EDES
        du =  [
            D[:GutGlucose](u,p,t,inputs)
            D[:PlasmaGlucose](u,p,inputs)
            D[:GlucoseIntegrator](u, p)
            D[:PlasmaInsulin](u, p, du[2])
            D[:RemoteInsulin](u, p)
            D[:NEFAInsulinDelays](u, p)
            D[:GutTriglyceride](u, p, t, inputs)
            D[:PlasmaTriglyceride](u, p, inputs)
            D[:PlasmaNEFA](u, p)]
        nothing
    end 
    sys!
end

simulation_specs = SimulationSpecs(
    initialvalues = (p, t0) -> [0., p.Gb, 0., p.Ib, 0., p.Ib, p.Ib, p.Ib, 0., 0., 0., p.TGb, p.NEFAb],
    timespan = (0., 480.)
)

# Define System
system_specs = SystemSpecs(
    merge(Dict(
        :NEFAInsulinDelays  => dNEFAInsulinDelays,
        :GutTriglyceride    => dTGgut,
        :PlasmaTriglyceride => dTGpl,
        :PlasmaNEFA         => dNEFApl,
    ),edes.architecture.equations),
    mealmodel_signature
)

MixedMealModel() = SysBioModel(
    system_specs,
    simulation_specs,
    parameters
)