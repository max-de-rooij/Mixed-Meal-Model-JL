using Distributed; addprocs(5, exeflags="--project")
using DelimitedFiles

# Base functionality
@everywhere begin
    using Plots, Measures, Statistics, StableRNGs

    # ODE Solver
    using DifferentialEquations

    # Benchmarking
    using BenchmarkTools

    include("../MixedMealModel.jl")
end

@everywhere function fit_individual(individual)
    rng = StableRNG(1847)
    glc, ins, trg, nfa, bwg, time = SampleData()
    model, constants, subject = ModelFromData(glc[individual,:], ins[individual,:], trg[individual,:], nfa[individual,:], bwg[individual])
    loss, simulationfunc = LossFunction(model, glc[individual,:], ins[individual,:], trg[individual,:], nfa[individual,:], Int.(time), constants, subject)
    #p_init = [1.35e-2, 3.80e-3, 5.82e-1, 0.00045, 0.0713, 208.88, 0.0163, 0.0119]
    p_lb = [0.005, 0, 0, 0, 0, 60., 0.005, 0]
    p_ub = [0.05, 1., 10., 1., 1., 720., 0.1, 1.]
    println("Fitting individual $(individual)")
    #println("Initial loss for individual $(individual): $(loss(p_init))")
    res, obj = FitModelLHC(loss, 5, p_lb, p_ub; rng = rng)
    return res, obj
end

@everywhere function FitAllIndividuals()
    @sync @distributed for i=1:5
        fit_individual(i)
    end
end

BenchmarkResult = @benchmark FitAllIndividuals() samples=100 seconds=60*60*2

# open("Windows_ShitMachine_Benchmark_ParamEstim.txt", "w+") do io
#     writedlm(io, BenchmarkResult.times)
# end

# BenchmarkResult