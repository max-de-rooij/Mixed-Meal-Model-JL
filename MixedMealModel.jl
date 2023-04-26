using DifferentialEquations

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

