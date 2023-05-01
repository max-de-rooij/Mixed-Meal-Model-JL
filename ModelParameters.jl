"""
Function to create a model parameter vector from the test subject. The parameters are by default fixed to their standard values, but can
be modified using their keyword arguments.

Arguments:
    sample_person   The test subject. A ComponentArray containing fields fasting_glucose, fasting_insulin, fasting_triglyceride, fasting_NEFA, body_weight, meal_glucose and meal_triglyceride

Optional arguments are the parameters. Read readme.md for an explanation of all parameters.
"""
function InitialParameters(person::NamedTuple; k1 = 0.0105, k2 = 0.28, k3 = 6.07e-3,
        k4 = 2.35e-4,
        k5 = 0.0424,
        k6 = 2.2975,
        k7 = 1.15,
        k8 = 7.27,
        k9 = 3.83e-2,
        k10 = 2.84e-1,
        sigma = 1.4,
        Km = 13.2,
        G_b = person.fasting_glucose,
        I_pl_b = person.fasting_insulin,
        G_liv_b = 0.043,
        spill = 30, 
        k11 = 0.00045, 
        ATL_max = 0.215, 
        K_ATL = 0.2, 
        k12 = 0.0713, 
        tau_LPL = 208.88, 
        k13 = 0.0088, 
        k14 = 0.0163, 
        k15 = 1e-5, 
        k16 = 0.0119
    )
    parameters = [
        k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, sigma, Km, G_b, I_pl_b, G_liv_b, spill, k11, ATL_max, K_ATL, k12, tau_LPL, k13, k14, k15, k16
    ]

    parameters
end

"""
Function to create a model initial_value vector from the test subject. The parameters are by default fixed to their standard values, but can
be modified using their keyword arguments.

Arguments:
    person   The test subject. A NamedTuple containing fields fasting_glucose, fasting_insulin, fasting_triglyceride, fasting_NEFA, body_weight, meal_glucose and meal_triglyceride

Optional arguments are the initial values. Read readme.md for an explanation of all initial values.
"""
function InitialValues(
        person::NamedTuple;
        mG_gut = 0., 
        G_pl = person.fasting_glucose, 
        G_int = 0., 
        I_pl = person.fasting_insulin,
        I_d₁ = 0.,
        I_d₂ = person.fasting_insulin,
        I_d₃ = person.fasting_insulin,
        I_d₄ = person.fasting_insulin, 
        NEFA_pl = person.fasting_NEFA, 
        mTG_gut₁ = 0., 
        mTG_gut₂ = 0.,
        mTG_gut₃ = 0.,
        TG_pl = person.fasting_triglyceride
    )
    
    [mG_gut, G_pl, G_int, I_pl, I_d₁, I_d₂, I_d₃, I_d₄, NEFA_pl, mTG_gut₁, mTG_gut₂, mTG_gut₃, TG_pl]
end

"""
Function to create a model constants vector from the parameter vector. The constants are by default fixed to their standard values, but can
be modified using their keyword arguments.

Arguments:
    parameters  a complete parameter vector made using the make_parameter_vector function.
    person   The test subject. A NamedTuple containing fields fasting_glucose, fasting_insulin, fasting_triglyceride, fasting_NEFA, body_weight, meal_glucose and meal_triglyceride

Optional arguments are the specific constants. Read readme.md for an explanation of all constants.
"""
function Constants(parameters, person;
        f_G = 0.005551,
        f_TG = 0.00113,
        f_I = 1.,
        V_G = (260/sqrt(person.body_weight/70))/1000,
        V_TG = (70/sqrt(person.body_weight/70))/1000,
        G_liv_b = parameters[15],
        tau_i = 31.,
        tau_d = 3.,
        G_threshold_pl = 9.,
        t_int = 30.,
        c1 = 0.1,
        c2 = parameters[15].*(parameters[12] + parameters[13])./parameters[13] - parameters[5].*parameters[15],
        c3 = parameters[7].*parameters[13]./(31 .* parameters[14]).*30
    )

    (f_G = f_G, f_TG = f_TG, f_I = f_I, V_G = V_G, V_TG = V_TG, G_liv_b = G_liv_b, τⁱ = tau_i, τᵈ = tau_d, G_threshold_pl = G_threshold_pl,
        t_int = t_int, c₁ = c1, c₂ = c2, c₃ = c3)
end

