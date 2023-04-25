function ModelLoss(model, inputs, glucose, insulin, triglyceride, nefa, timepoints)

    problem = ODEProblem(model, inputs)
    fixed_parameters = non_estimated_params(model)

    timepoints = Vector{Int}(timepoints)

    function Loss(p)

        # Solve ODE
        pp = [p; fixed_parameters]
        s = DifferentialEquations.solve(problem, KenCarp4(), p=pp, saveat=1, sensealg=ForwardDiffSensitivity(), verbose=false)
        sol = Array(s)

        # glucose error component
        glucose_error = sol[2,timepoints.+1].-glucose
        glucose_loss = sum(abs2, glucose_error)/maximum(glucose)

        # insulin error component
        insulin_error = sol[4,timepoints.+1].-insulin
        insulin_loss = sum(abs2, insulin_error)/maximum(insulin)

        # TG error component
        TG_error = sol[12,timepoints.+1].-triglyceride
        TG_loss = sum(abs2, TG_error)/maximum(triglyceride)

        # NEFA error component
        NEFA_error = sol[13,timepoints.+1].-nefa
        NEFA_loss = sum(abs2, NEFA_error)/maximum(nefa)

        # Fit error
        scaling_term = maximum(glucose)
        fit_error = scaling_term.*(glucose_loss + insulin_loss + TG_loss + NEFA_loss)

        # error if AUC of roc of gut glucose < meal content
        G_gut = t -> (pp.k2 * (pp.f / (pp.vg * inputs[2]))) .* sol[1, t .+ 1]
        AUC_G = (sum(G_gut, 1:239) + (G_gut(0) + G_gut(240))/2) * ((pp.vg * inputs[2])/pp.f)
        err_AUC_G = abs(AUC_G - inputs[1])/10000

        # error if AUC of roc of TG in plasma < meal content
        TG_gut = t -> (pp.k14 * (pp.fTG/(pp.vTG*inputs[2]))) .* sol[11, t .+ 1]
        AUC_TG = (sum(TG_gut, 1:479) + (TG_gut(0) + TG_gut(480))/2) * ((pp.vTG*inputs[2])/pp.fTG)
        err_AUC_TG = abs(AUC_TG - inputs[3])/10000
    
        # constrain steady state TG to measured fasting value
        TG_steady_state = (pp.TGb - sol[11, 721])[1]

        # constrain steady state NEFA to measured fasting value
        model_fasting_NEFA = (3*(pp.fspill/100)*pp.k11*pp.TGb*pp.Ib + pp.ATLmax/(1+pp.KATL*(pp.Ib)^2)) / pp.k12
        NEFA_diff = pp.NEFAb - model_fasting_NEFA

        # Regularisation error
        regularisation_error = err_AUC_G + err_AUC_TG + abs2(TG_steady_state) + abs2(NEFA_diff)

        # Combined Loss Value
        fit_error + regularisation_error
    end

    function OutputValues(p)
        # Solve ODE
        pp = [p; fixed_parameters]
        DifferentialEquations.solve(problem, Tsit5(), p=pp, saveat=1)
    end
        

    Loss, OutputValues
end

function FitModelLHC(optf, n, lb, ub, ϵ = 1e-9)
    parameter_sets = LHCoptim(n, length(lb), 1000)[1] ./ n
    # scale parameter sets
    ubx = (1-ϵ).*ub;
    lbx = (1+ϵ).*lb;
    parameter_sets = (parameter_sets' .* (ubx.-lbx)) .+ lbx

    results = Vector{Float64}[]
    objectives = Float64[]

    callback = (p, l) -> begin 
        println("Loss: $(l)")
        false
    end
    
    for it in axes(parameter_sets,2)
        optprob = OptimizationProblem(optf, parameter_sets[:, it],lb=lb, ub=ub)
        local_optimizer = Optim.LBFGS(linesearch=BackTracking(order=3))
        try
            sol = Optimization.solve(optprob, local_optimizer, x_tol=1e-8, f_tol = 1e-6, g_tol=1e-6, f_calls_limit=1000)
            push!(results, sol.minimizer)
            push!(objectives, sol.objective)
        catch 
            println("Optimization Failed... Resampling...")
        end
    end



    return results, objectives
end