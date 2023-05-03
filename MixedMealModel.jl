using DifferentialEquations, SciMLSensitivity, LatinHypercubeSampling, Optimization, OptimizationOptimJL, LineSearches
include("ModelParameters.jl")

"""
Obtain the ODE function based on the constants and the meal input data specified. The resulting function can be used to
create an ODEProblem object using the DifferentialEquations.jl package.

Arguments:
    constants   contains the model constants, can be easily defined with the make_constants function after the initial parameter vector is known
    input_data  contains the test subject specific data, such as the body weight and the steady-state plasma levels of glucose, insulin, TG and NEFA.
    
The output is a function with the signature f(D, u, p, t), where
    D   derivatives
    u   state variables
    p   parameters
    t   time
"""
function MixedMealModel(constants::NamedTuple, input_data::NamedTuple)

    system! = function (D, u, p, t)
        # define model state variables
        mG_gut, G_pl, G_int, I_pl, I_d₁, I_d₂, I_d₃, I_d₄, NEFA_pl, mTG_gut₁, mTG_gut₂, mTG_gut₃, TG_pl = u

        # model input
        D_meal_G = input_data.meal_glucose
        D_meal_TG = input_data.meal_triglyceride
        BW = input_data.body_weight

        # model parameters
        k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, σ, Kₘ, G_b, I_pl_b, G_liv_b, spill, k11, ATL_max, K_ATL, k12, τ_LPL, k13, k14, k15, k16 = p
        
        # model equations

        # appearance of glucose from the meal 
        G_meal = σ*(k1.^σ)*t.^(σ-1) .* exp(-1*(k1.*t).^σ) * D_meal_G
        D[1] = G_meal - k2*mG_gut

        # plasma glucose equations
        # glucose flux accross the liver
        G_liv = G_liv_b - k4.*constants.f_I.*I_d₁ - k3.*(G_pl-G_b)
        # glucose mass in gut
        G_gut = k2.*(constants.f_G/(constants.V_G*BW)).*mG_gut
        # insulin-independent utilization of glucose
        U_ii = G_liv_b*((Kₘ + G_b)./G_b).*(G_pl./(Kₘ + G_pl))
        # insulin-dependent utilization of glucose
        U_id = k5.*constants.f_I.*I_d₁.*(G_pl./(Kₘ + G_pl))
        # renal extraction of plasma glucose
        U_ren = (constants.c₁./(constants.V_G*BW).*(G_pl - constants.G_threshold_pl))*(G_pl > constants.G_threshold_pl)

        D[2] = G_liv + G_gut - U_ii - U_id - U_ren

        # plasma insulin equations
        # integrator term of the PID
        D[3] = G_pl - G_b
        # pancreatic production of insulin
        I_pnc = (constants.f_I.^-1).*(k6.*(G_pl - G_b) + (k7/constants.τⁱ).*G_int + (k7/constants.τⁱ).*G_b + (k8.*constants.τᵈ).*D[2])
        # liver insulin
        I_liv = k7.*(G_b./(constants.f_I.*constants.τⁱ.*I_pl_b)).*I_pl
        # transport of insulin from plasma to remote
        I_rem = k9*(I_pl-I_pl_b)
        # plasma insulin roc
        D[4] = I_pnc - I_liv - I_rem
        # remote insulin
        D[5] = k9*(I_pl - I_pl_b) - k10*I_d₁

        # insulin delays for NEFA
        D[6] = (3/τ_LPL).*(I_pl - I_d₂)
        D[7] = (3/τ_LPL).*(I_d₂ - I_d₃)
        D[8] = (3/τ_LPL).*(I_d₃ - I_d₄)

        # uptake of circulating TG into tissues through LPL lipolysis
        LPL = k11.*TG_pl.*I_d₄

        # NEFA
        frac_spill = (1/100).*(spill.*(I_pl_b./I_d₂))
        # plasma NEFA
        D[9] = 3 .* frac_spill.*LPL + ATL_max./(1+K_ATL.*(I_d₂).^2) - k12.*NEFA_pl;

        # appearance of TG from meal
        TG_meal = σ*(k13.^σ)*t.^(σ-1) .* exp(-1*(k13.*t).^σ) * D_meal_TG
        # TG_mass in gut
        D[10] = TG_meal - k14*mTG_gut₁
        D[11] = k14*mTG_gut₁ - k14*mTG_gut₂
        D[12] = k14*mTG_gut₂ - k14*mTG_gut₃

        # plasma TG
        # TG from meal (chylomicron)
        TG_gut = k14.*(constants.f_TG/(constants.V_TG*BW)).*mTG_gut₃
        
        # secretion of endogenous TG from the liver (VLDL)
        VLDL = k16-k15.*(I_d₄-I_pl_b)
        # roc of plasma TG
        D[13] = VLDL + TG_gut - LPL
        # steady-state: k16 = k11.*TG_pl_b.*I_pl_b 

    end
    return system!
end

"""
Get the output of the model from the solution of the ODE System.

Arguments:
    solution        ODE solution from the solve function
    parameters      Model parameter vector
    constants       Model constants
    test_subject    Model test subject
    time (optional) Optional time vector different than the solution
    
The output is a NamedTuple containing:
    - Gut glucose
    - Hepatic glucose flux
    - Glucose uptake by tissue
    - Plasma glucose
    - Plasma insulin
    - Gut triglyceride
    - Liver triglyceride
    - Plasma triglyceride
    - Plasma NEFA
"""
function ModelOutput(solution, parameters, constants::NamedTuple, test_subject::NamedTuple; time = solution.t)

        # glucose in the gut
        mG_gut = solution(time, idxs=1) #[k.mG_gut for k in solution.u]
        G_gut = parameters[2].*(constants.f_G/(constants.V_G*test_subject.body_weight)).*mG_gut
        
    
        # hepatic glucose flux
        k3 = parameters[3]
        k4 = parameters[4]
        G_liv_b = parameters[15]
        G_b = parameters[13]
        I_d₁ = solution(time, idxs=5) #[k.I_d₁ for k in solution.u]
        G_pl = solution(time, idxs=2) #[k.G_pl for k in solution.u]
        G_liv = G_liv_b .- k4.*constants.f_I.*I_d₁ .- k3.*(G_pl.-G_b)
    
        # glucose uptake into tissue
        KM = parameters[12];
        k5=parameters[5];
        U_ii = G_liv_b*((KM .+ G_b)./G_b).*(G_pl./(KM .+ G_pl));
    
        U_id = k5.*constants.f_I.*I_d₁.*(G_pl./(KM .+ G_pl));
    
        G_U = U_ii .+ U_id;
    
    
        # Plasma glucose
    
        # plasma insulin
        I_pl = solution(time, idxs=4) # [k.I_pl for k in solution.u]
    
        # TG from gut
        k14 = parameters[23]
        M_TG_gut = solution(time, idxs=12) #[k.mTG_gut₃ for k in solution.u]
        TG_gut = k14.*(constants.f_TG./(constants.V_TG.*test_subject.body_weight)).*M_TG_gut
    
        # TG from liver
        I_d₄ = solution(time, idxs=8) #[k.I_d₄ for k in solution.u]
        I_pl_b = parameters[14]
        k15 = parameters[24]
        k16 = parameters[25]
    
        VLDL = k16 .- k15.*(I_d₄.-I_pl_b)
    
        # plasma TG
        TG_pl = solution(time, idxs=13) #[k.TG_pl for k in solution.u]
    
        # Plasma NEFA
        NEFA_pl = solution(time, idxs=9) #[k.NEFA_pl for k in solution.u]
    
        return (
            gut_glucose = G_gut,
            hepatic_glucose_flux = G_liv,
            glucose_tissue_uptake = G_U,
            plasma_glucose = G_pl,
            plasma_insulin = I_pl,
            gut_TG = TG_gut,
            liver_TG = VLDL,
            plasma_TG = TG_pl,
            plasma_NEFA = NEFA_pl
        )
    end


"""
Obtain sample data for five subjects for fitting the Mixed Meal Model

Arguments:
    -
    
Outputs:
    - Glucose (5x8)
    - Insulin (5x8)
    - Triglyceride (5x8)
    - NEFA (5x8)
    - Body Weight (5x1)
    - Timepoints (8x1)
"""
function SampleData()

    glc = [4.97828626648081	9.83226219754181	6.91736265000290	3.99736583151076	4.74531292225800	5.43827322143290	5.00061029702766	4.90628637452010;
    4.44438211792684	8.56497936784908	5.63459914693781	2.87846029141051	4.38475106119756	4.81819939162605	4.68440337786608	5.03050899831009;
    5.86649175256700	9.76226946026824	6.96041362059030	5.75415042357602	5.49442472724051	6.04402280406414	6.20665742418169	6.05245308331924;
    4.05422828334834	6.68317282007679	5.51316702906051	4.38821628507995	4.08278484538879	4.32903985847918	3.87347487223252	4.14703615153245;
    5.09861020045252	7.19169651136166	6.82073905407388	4.41846994550959	4.35813489280548	5.26064945251012	5.11261281908290	4.95616017748044
    ]
    ins = [21.7542994608033	287.361969174996	294.323541530305	85.4388252908261	18.0684917209633	21.0182483747640	23.1679411923448	21.8618244436204;
    12.5404266739517	124.050624391434	112.050042225095	38.5023639157180	4.04338539971474	13.9239709608373	13.4555353291966	14.3638940482423;
    14.2802823649251	187.522140183945	134.487839622436	29.1762282256049	14.4252356798284	13.8846272694087	14.1392062008499	14.4034973030377;
    19.8623112280412	199.253276210798	258.747486618497	108.466775252880	41.2783085983230	28.2728482488501	22.5314861118322	21.0737318920216;
    108.760022278870	163.060430670949	207.080258023449	185.748135553145	125.786642682319	109.125617637325	104.745417655491	97.2648413826163  
    ]
    trg = [1.16871366007653	1.58440178335270	1.78195011372397	2.39301829368748	2.41622964225137	2.31174933441785	1.93219870120156	1.72103932334563;
    1.12064493587155	0.988185237587493	1.12903104782257	1.27025369534736	1.38325830618782	1.70424182674534	1.74785291292704	1.58034015443830;
    0.934147873129793	0.948495684902739	1.06182083372904	1.19875606178644	1.15075004830806	0.964931259331579	0.959197106220412	0.941575651575427;
    0.911568878711074	0.917659723900865	1.04804850665628	1.06069565177322	1.33975089951343	1.32507244917361	1.25922207558199	0.954325323451108;
    1.25060974085694	1.17562009272981	1.09605143482126	0.948511358144563	1.12086592905305	1.38356985274728	1.74007908698886	1.96096980463590
    ]
    nfa = [0.335735045292724	0.306808022891325	0.175773322340940	0.126453528770619	0.186422309973495	0.317877874980899	0.504623299086489	0.561199670138616;
    0.497614523572355	0.214246896777152	0.0965385883289903	0.0679221665148035	0.210855464643708	0.457945252128255	0.576852516292241	0.515581280746099;
    0.156640771311381	0.0771447415470159	0.0320943309717143	0.0831605186940559	0.189924558213506	0.296184914024023	0.296916721601158	0.215514420118521;
    0.499717860234900	0.311498736265132	0.175340122607563	0.108305580613699	0.127238842499880	0.206930006398300	0.340526704496308	0.388220201643060;
    0.323799139072395	0.297458153444511	0.252292171822141	0.192210379268320	0.192298067008140	0.314571176721740	0.495631829890529	0.471625083527584
    ]

    bwg = [85, 76, 71, 91, 120]
    time = [0.,30.,60.,120.,180.,240.,360.,480.]

    return glc, ins, trg, nfa, bwg, time
end


"""
Create an ODEProblem from the input data of one subject.

Arguments:
    - Glucose: Glucose data
    - Insulin: Insulin data 
    - Triglyceride: TG data
    - NEFA: NEFA data
    - BodyWeight: Body Weight
    Optional:
    - TimeSpan: simulation timespan, default [0, 720]
    - MealG: Meal Glucose, default 75000. mg
    - MealTG: Meal TG, default 60000. mg

Output:
    - ODEProblem containing the model and the default parameters

"""
function ModelFromData(Glucose, Insulin, Triglyceride, NEFA, BodyWeight; TimeSpan = (0., 720.), MealG = 75000., MealTG = 60000.)

    sample_person = (
        fasting_glucose = Glucose[1],
        fasting_insulin = Insulin[1],
        fasting_triglyceride = Triglyceride[1],
        fasting_NEFA = NEFA[1],
        body_weight = BodyWeight,
        meal_glucose = MealG,
        meal_triglyceride = MealTG
    )

    parameters = InitialParameters(sample_person, k1 = 0.0164, k5=0.0564, 
    k6 = 2.7341, k11 = 0.00035, k12 = 0.0822, tau_LPL = 187.88, 
    k14 = 0.0392, k16 = 0.0135)

    constants = Constants(parameters, sample_person);

    tspan = TimeSpan
    u0 = InitialValues(sample_person)
    MealModel! = MixedMealModel(constants, sample_person)
    return ODEProblem(MealModel!, u0, tspan, parameters), constants, sample_person
end

function LossFunction(Model, Glucose, Insulin, Triglyceride, NEFA, TimePoints, Constants, Subject)

    initialparameters = Model.p;

    function Loss(p)
        
        parametervector = [
            p[1]; 
            initialparameters[2:4]; 
            p[2:3]; 
            initialparameters[7:16]; 
            p[4]; 
            initialparameters[18:19]; 
            p[5:6]; 
            initialparameters[22]; 
            p[7]; 
            initialparameters[24]; 
            p[8]
        ]

        s = solve(Model, p=parametervector, saveat=[0:480; 720], sensealg=ForwardDiffSensitivity(), verbose=false)
        sol = Array(s)

        size(sol, 2) == 482 || throw(ErrorException("ODE Solver Timestep is too small!"))

        # glucose error component
        glucose_error = (sol[2,TimePoints.+1].-Glucose)./maximum(Glucose)

        # insulin error component
        insulin_error = (sol[4,TimePoints.+1].-Insulin)./maximum(Insulin)

        # TG error component
        TG_error = (sol[13,TimePoints.+1].-Triglyceride)./maximum(Triglyceride)

        # NEFA error component
        NEFA_error = (sol[9,TimePoints.+1].-NEFA)./maximum(NEFA)

        # Fit error
        scaling_term = maximum(Glucose)
        fit_error = scaling_term.*[glucose_error; insulin_error; TG_error; NEFA_error]

        # error if AUC of roc of gut glucose < meal content
        G_gut = t -> (parametervector[2] .* (Constants.f_G ./ (Constants.V_G .* Subject.body_weight))) .* sol[1, t .+ 1]
        AUC_G = (sum(G_gut, 1:239) + (G_gut(0) + G_gut(240))/2) * ((Constants.V_G * Subject.body_weight)/Constants.f_G)
        err_AUC_G = abs(AUC_G - Subject.meal_glucose)/10000

        # error if AUC of roc of TG in plasma < meal content
        TG_gut = t -> (parametervector[23] .* (Constants.f_TG./(Constants.V_TG.*Subject.body_weight))) .* sol[12, t .+ 1]
        AUC_TG = (sum(TG_gut, 1:479) + (TG_gut(0) + TG_gut(480))/2) * ((Constants.V_TG*Subject.body_weight)/Constants.f_TG)
        err_AUC_TG = abs(AUC_TG - Subject.meal_triglyceride)/10000

        # constrain steady state glucose to measured fasting value
        G_steady_state = (Subject.fasting_glucose - sol[2, 301])[1]
    
        # constrain steady state TG to measured fasting value
        TG_steady_state = (Subject.fasting_triglyceride - sol[13, 482])[1]

        # constrain steady state NEFA to measured fasting value
        model_fasting_NEFA = (3*(parametervector[16]/100)*parametervector[17]*Subject.fasting_triglyceride*Subject.fasting_insulin + parametervector[18]/(1+parametervector[19]*(Subject.fasting_insulin)^2)) / parametervector[20]
        NEFA_diff = Subject.fasting_NEFA - model_fasting_NEFA

        # non-negative VLDL flux
        VLDL_nonneg = sum(abs2, min.(0, parametervector[25] .- parametervector[24].*(sol[8,:].-Subject.fasting_insulin)))

        # Regularisation error
        regularisation_error = [err_AUC_G, err_AUC_TG, G_steady_state, VLDL_nonneg, TG_steady_state, 8 .*NEFA_diff]

        # Combined Loss Value
        sum(abs2, [fit_error; regularisation_error])
    end

    function SimulateModel(p)
        parametervector = [
            p[1]; 
            initialparameters[2:4]; 
            p[2:3]; 
            initialparameters[7:16]; 
            p[4]; 
            initialparameters[18:19]; 
            p[5:6]; 
            initialparameters[22]; 
            p[7]; 
            initialparameters[24]; 
            p[8]
        ]

        solve(Model, p=parametervector)
    end

    return Loss, SimulateModel
end


function FitModelLHC(loss, n, lb, ub; ϵ = 1e-9, rng = StableRNG(1234))

    optf = OptimizationFunction((x,λ) -> loss(x), Optimization.AutoZygote())
    parameter_sets = LHCoptim(n, length(lb), 1000; rng = rng)[1] ./ n
    # scale parameter sets
    ubx = (1-ϵ).*ub;
    lbx = (1+ϵ).*lb;
    parameter_sets = (parameter_sets' .* (ubx.-lbx)) .+ lbx

    results = []
    objectives = Float64[]

    callback = (p, l) -> begin 
        println("Loss: $(l)")
        false
    end
    
    for it in axes(parameter_sets,2)
        optprob = OptimizationProblem(optf, parameter_sets[:, it],lb=lb, ub=ub)
        local_optimizer = Optim.LBFGS(linesearch=BackTracking(order=3))
        try
            sol = Optimization.solve(optprob, local_optimizer, f_calls_limit=1000, x_tol=1e-8, f_tol = 1e-6, g_tol=1e-6, maxiters=400)
            push!(results, sol.minimizer)
            push!(objectives, sol.objective)
        catch 
            println("Optimization Failed... Resampling...")
        end
    end

    return results, objectives
end


    
