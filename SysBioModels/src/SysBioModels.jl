module SysBioModels

    using ComponentArrays
    using SciMLBase
    import DiffEqBase.ODEProblem
    using OrdinaryDiffEq
    import Base.+

    mutable struct SimulationSpecs{T<:Real, U<:Union{Nothing, OrdinaryDiffEqAlgorithm}}
        initialvalues::Union{Function, <:AbstractVector{T}}
        timespan::Tuple{T, T}
        alg::U
        function SimulationSpecs(;initialvalues::Union{Function, <:AbstractVector{T}}, timespan::Tuple{T, T}, algorithm::U = nothing) where {T<:Real, U<:Union{Nothing, OrdinaryDiffEqAlgorithm}}
            return new{T, U}(initialvalues, timespan, algorithm)
        end
    end

    function DefaultSignature(Equations::Dict{Symbol, Function})
        sys! = function(du, u, p, t)
            for i in eachindex(values(Equations))
                du[i] = values(Equations)[i](du, u, p, t)
            end
        end
        sys!
        end

    mutable struct SystemSpecs
        equations::Dict{Symbol, Function}
        systemsignature::Function
        function SystemSpecs(equations::Dict{Symbol, Function}, systemsignature::Function=DefaultSignature)
            return new(equations, systemsignature)
        end
    end

    mutable struct ModelParameters{T}
        estimated::ComponentVector{T}
        fixed::ComponentVector{T}
        constants::ComponentVector{T}
        function ModelParameters(;estimated::ComponentVector{T} = ComponentVector(), fixed::ComponentVector{T} = ComponentVector(), constants::ComponentVector{T} = ComponentVector()) where T
            return new{T}(estimated, fixed, constants)
        end
    end

    struct SysBioModel{T, U}
        architecture::SystemSpecs
        simulator::SimulationSpecs{T, U}
        parameters::ModelParameters{T}
    end

    function system(model::SysBioModel, args...)
        model.architecture.systemsignature(model.architecture.equations, args...)
    end

    function params(model::SysBioModel)
        [[model.parameters.fixed; model.parameters.estimated]; model.parameters.constants]
    end

    function estimated_params(model::SysBioModel)
        model.parameters.estimated
    end

    function non_estimated_params(model::SysBioModel)
        [model.parameters.fixed; model.parameters.constants]
    end

    function initials(model::SysBioModel)
        model.simulator.initialvalues
    end

    function timespan(model::SysBioModel)
        model.simulator.timespan
    end

    function ODEProblem(model::SysBioModel, args...)
        ODEProblem(system(model, args...), initials(model), timespan(model), params(model))
    end

    function +(a::ModelParameters{T}, b::ModelParameters{T}) where T
        ModelParameters(
            estimated = [a.estimated; b.estimated],
            fixed     = [a.fixed; b.fixed],
            constants = [a.constants; b.constants]
        )
    end

    function simulate(model::SysBioModel{T, Nothing}, args...) where T
        solve(ODEProblem(model, args...), save_everystep=false)
    end

    function simulate(model::SysBioModel{T, U}, args...) where {T, U<:OrdinaryDiffEqAlgorithm}
        solve(ODEProblem(model, args...), model.simulator.alg, save_everystep=false)
    end

    # Custom structs
    export SimulationSpecs, SystemSpecs, ModelParameters, SysBioModel
    # Custom helper functions
    export system, params, estimated_params, non_estimated_params, initials, timespan 
    # Implementations of other packages
    export ODEProblem, +
    # Workflow
    export simulate
end