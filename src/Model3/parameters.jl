
"""
Global parameters with defaults.
"""
glb() = ComponentVector(
    t_max = 3 * 365., # max simulation time (real simulation time can be shorter when callbacks are triggered) [d]
    N0 = 1., # initial number of individuals (always 0 for pure-ODE system) [#]
    food_dynamic = 1., # boolean indicator of whether food abundance is treated as a dynamic variable [-]
    dX_in_aq = 20., # food input rate in aquatic medium [mg d^-1]
    dX_in_ter = 20., # food input rate in terrestric environment [mg d^-1]
    k_V_aq = 0., # resource dilution rate in the aquatic medium [d^-1] 
    k_V_ter = 0., # resource dilution rate in the terrestric medium (interpretable as background mortality of the prey) [d^-1]
    V_patch_aq = 1., # simulated volume of the aquatic environment [L] 
    A_patch_ter = 1., # simulated area of the terrestric environment [m^2]
    T_aq = 293.15, # ambient temperature in aquatic environment [K]
    T_ter = 293.15, # ambient temperature in terrestric environment [K]
    C_W1 = 0., # aquatic exposure concentrations of substance 1 [e.g. μg/L]
    C_W2 = 0., # aquatic exposure concentrations of substance 2 [e.g. μg/L]
    pathogen_inoculation_dose = 0., # amount of pathogen spores added to aquatic medium [# spores]
    pathogen_inoculation_time = 30., # time-point of pathogen inoculation [d]
)

"""
Species-specific parameters with defaults.
"""
spc() = ComponentVector(

    #=
    Metaparameters
    =#

    Z = Dirac(1.), # mass-based zoom factor [-]
       
    #=
    Physiological baseline (DEB) parameters
    =#

    X_emb_int = 1, # Initial vitellus (≈ dry mass of an egg) [mg]
    K_X_lrv = 1.,  # larval half-saturation constant for food uptake [mg L^-1]
    K_X_juv = 1., # juvenile and adult half-saturation constant for food uptake [mg m^-2]
    dI_max_emb = 1, # embryonic maximum specific ingestion rate [mg mg^-2/3 d^-1]
    dI_max_lrv = 1, # larval maximum specific ingestion rate [mg mg^-2/3 d^-1]
    dI_max_juv = 1, # juvenile and adult maximum specific ingestion rate [mg mg^-2/3 d^-1]
    kappa_emb = 0.8, # embryonic to metamorph allocation fraction; default value is suggested value from DEBkiss book []
    kappa_juv = 0.8, # juvenile and adult allocation fraction to soma; default value is suggested value from DEBkiss book [-]
    eta_IA = 0.54, # assimilation efficiency; default value from Pfab et al. (2020) [-]
    eta_AS_emb = 0.4, # embryonic to metamorph growth efficiency; default value from Pfab et al. (2020) [-]
    eta_AS_juv = 0.4, # juvenile and adult growth efficiency [-]
    eta_AR = 0.95, # reproduction efficiency; default value is suggested value from DEBkiss book [-]
    eta_SA = 0.8, # shrinking efficiency [-]
    k_M_emb = 0.11, # embryonic to metamorph somatic maintenance rate constant; default value from Pfab et al. (2020)
    k_M_juv = 0.11, # juvenile and adult somatic maintenance rate constant; initially assumed equal to k_M_emb
    k_J_emb = 0.027, # embryonic to metamorph maturity maintenance rate constant; initially calculated based on (1-kappa)/kappa-ratio [d^-1]
    k_J_juv = 0.027, # juvenile and adult maturity maintenance rate constant; initially assumed equal to k_J_emb [d^-1]
    k_C = 0.1, # metamorph mobilization rate of reserve buffer [d^-1]
    H_j1 = 1, # maturity (resreve buffer) at the start of metamorphosis (decline of feeding) [mg]
    H_j2 = 1e-3, # maturity at the end of metamorphosis [mg]
    H_p = 55., # maturity at puberty [mg]

    T_A = 8000., # Arrhenius temperature [K]
    T_ref = 293.15, # reference 
    b_T = 40., # effect strength of temperature on resource allocation

    h_b_aq = 1e-3, # background mortality aquatic stages incl. metamorphs [d^-1]
    h_b_ter = 1e-3, # background mortality terrestric stages [d^-1]

    watercontent_larvae = 0.9, # water content of larvae [-]
    watercontent_juveniles = 0.75, # water content of juveniles and adults [-]

    aux = ComponentVector( # axiliary parameters
        tau_R = 365., # reproduction period [d]
        h_S = 1e-3, # starvation hazard rate [d^-1]
        S_rel_crit = 1/3, # critical fraction of structure that can be lost before h_S is triggered [-]
        a_max = Truncated(Normal(3650, 365), 0, Inf) # maximum life span [d]
    ), 
    tktd = tktdparams_binary_loglogistic() # tktd parameters 
)

function tktdparams_binary_loglogistic()
    return ComponentVector(
        k_D1_G = 0.,
        k_D1_M = 0.,
        k_D1_A = 0.,
        k_D1_R = 0.,
        k_D1_H = 0., # increase in H_j1
        k_D1_κ_pos = 0., # increase in κ
        k_D1_κ_neg = 0., # decrease in κ

        e1_G = 1.,
        e1_M = 1.,
        e1_A = 1.,
        e1_R = 1.,
        e1_H = 1.,
        e1_κ_pos = 1., # increase in κ
        e1_κ_neg = 1., # decrease in κ

        b1_G = 2.,
        b1_M = 2.,
        b1_A = 2.,
        b1_R = 2.,
        b1_H = 2., # increase in H_j1
        b1_κ_pos = 2., # increase in κ
        b1_κ_neg = 2., # decrease in κ

        k_D2_G = 0.,
        k_D2_M = 0.,
        k_D2_A = 0.,
        k_D2_R = 0.,
        k_D2_H = 0., # increase in H_j1
        k_D2_κ_pos = 0., # increase in κ
        k_D2_κ_neg = 0., # decrease in κ

        e2_G = 1.,
        e2_M = 1.,
        e2_A = 1.,
        e2_R = 1.,
        e2_H = 1., # increase in H_j1
        e2_κ_pos = 1., # increase in κ
        e2_κ_neg = 1., # decrease in κ

        b2_G = 2.,
        b2_M = 2.,
        b2_A = 2.,
        b2_R = 2.,
        b2_H = 2., # increase in H_j1
        b2_κ_pos = 2., # increase in κ
        b2_κ_neg = 2., # decrease in κ
    )
end


# we can use the same parameter vector for the linear model, 
# interpreting `e` as threshold and `b` as slope
const tktdparams_binary_linear = tktdparams_binary_loglogistic

defaultparams() = ComponentVector(glb = glb(), spc = spc())

function generate_individual_params(p; kwargs...)

    ind = EcotoxSystems.getval.(p.spc) |> 
    x -> Real.(x) |> 
    x -> begin
        x.dI_max_emb *= x.Z^(1/3)
        x.dI_max_lrv *= x.Z^(1/3)
        x.dI_max_juv *= x.Z^(1/3)
        x.H_p *= x.Z
        x.X_emb_int *= x.Z
        x.K_X_lrv *= x.Z^(1/3)
        x.K_X_juv *= x.Z^(1/3)
        x
    end

    return ComponentVector(
        glb = p.glb, 
        ind = ind; 
        kwargs...
    )

end