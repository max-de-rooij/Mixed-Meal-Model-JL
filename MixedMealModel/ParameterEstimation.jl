using Distributed; addprocs(5)
@everywhere begin
    import Pkg
    Pkg.activate(".")
end

@everywhere begin
    using DifferentialEquations, BenchmarkTools

    using Optimization, OptimizationOptimJL, OptimizationOptimisers, SciMLSensitivity, Zygote, LineSearches

    using Plots, SharedArrays

    using LatinHypercubeSampling

    using EdesVersions, SysBioModels, ComponentArrays
    include("MealModel.jl")
    include("BenchmarkingUtils.jl")
    include("ParamEstimUtils.jl")


    # Load the data
    glc, ins, trg, nfa, bwg, time = LoadSampleData()
end

function FitAllIndividuals()
    @sync @distributed for individual in eachindex(bwg)
        inputs = [
        75000., # Meal Glucose
        bwg[individual],    # Body Mass
        60000., # Meal TG
        ]

        # Load the base model
        model = MixedMealModel()

        model.simulator.timespan = (0., 720.)

        model.parameters.fixed.Gb = glc[individual, 1]
        model.parameters.fixed.Ib = ins[individual, 1]
        model.parameters.fixed.TGb = trg[individual, 1]
        model.parameters.fixed.NEFAb = nfa[individual, 1]
        p = non_estimated_params(model)
        model.simulator.initialvalues = [0., p.Gb, 0., p.Ib, 0., p.Ib, p.Ib, p.Ib, 0., 0., 0., p.TGb, p.NEFAb]
        loss, simulator = ModelLoss(model, inputs, glc[individual, :], ins[individual, :], trg[individual, :], nfa[individual, :], time)

        println("Initial loss for individual $(individual): $(loss(model.parameters.estimated))")

        optf = OptimizationFunction((x,λ) -> loss(x), Optimization.AutoZygote())

        lb = ComponentVector(k11 = 1e-6, k12 =1e-6, tauLPL = 60., k14 = 0.005, k16 = 1e-6, k1 = 0.005, k5 = 1e-6, k6 = 1e-6)
        ub = ComponentVector(k11 = 1., k12 = 1., tauLPL = 600., k14 = 0.1, k16 = 1., k1 = 0.1, k5 = 1., k6 = 5.)

        FitModelLHC(optf, 25, lb, ub)
    end
    nothing
end

benchmark = @benchmarkable FitAllIndividuals() samples=100 seconds=60*60*2
run(benchmark)
# println("Individual $(individual): $(benchmarkresult)")

# glc_plot = scatter(time, glc[individual, :], labels="Data", title="Plasma Glucose", ylabel="Concentration [mmol/L]")
# plot!(glc_plot, outs, idxs=2, labels="Model")
# ins_plot = scatter(time, ins[individual, :], labels="", title="Plasma Insulin", ylabel="Concentration [uIU/mL]")
# plot!(ins_plot, outs, idxs=4, labels="")
# trg_plot = scatter(time, trg[individual, :], labels="", title="Plasma TG", ylabel="Concentration [mmol/L]")
# plot!(trg_plot, outs, idxs=12, labels="")
# nfa_plot = scatter(time, nfa[individual, :], labels="", title="Plasma NEFA", ylabel="Concentration [mmol/L]")
# plot!(nfa_plot, outs, idxs=13, labels="")

# output_plot = plot(glc_plot, ins_plot, trg_plot, nfa_plot)