
"""
Global parameters with defaults.
"""
glb = ComponentVector(
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
spc = ComponentVector(

    #=
    Metaparameters
    =#

    Z = Dirac(1.), # zoom factor
       
    #=
    Physiological baseline (DEB) parameters
    =#

    X_emb_int = 1, # Initial vitellus (≈ dry mass of an egg)
    K_X_lrv = 1.,  # larval half-saturation constant for food uptake
    K_X_juv = 1., # juvenile and adult half-saturation constant for food uptake
    dI_max_emb = 1, # embryonic maximum specific ingestion rate
    dI_max_lrv = 1, # larval maximum specific ingestion rate
    dI_max_juv = 1, # juvenile and adult maximum specific ingestion rate
    kappa_emb = 0.8, # embryonic to metamorph allocation fraction; default value is suggested value from DEBkiss book
    kappa_juv = 0.8, # juvenile and adult allocation fraction to soma; default value is suggested value from DEBkiss book
    gamma = 0.5, # larval allocation fraction to metamorphic reserves (split occurs downstream of κ-split)
    eta_IA = 0.54, # assimilation efficiency; default value from Pfab et al. (2020)
    eta_AS_emb = 0.4, # embryonic to metamorph growth efficiency; default value from Pfab et al. (2020)
    eta_AS_juv = 0.4, # juvenile and adult growth efficiency; initially assumed eval to eta_AS_emb
    eta_AR = 0.95, # reproduction efficiency [-]; default value is suggested value from DEBkiss book and also the default in add-my-pet
    eta_SA = 0.8, # shrinking efficiency
    k_M_emb = 0.11, # embryonic to metamorph somatic maintenance rate constant; default value from Pfab et al. (2020)
    k_M_juv = 0.11, # juvenile and adult somatic maintenance rate constant; initially assumed equal to k_M_emb
    delta_k_M_mt = 1., # metamorph somatic maintenance rate constant, relative to that of embryos andd larvae (not used)
    k_J_emb = 0.027, # embryonic to metamorph maturity maintenance rate constant; initially calculated based on (1-kappa)/kappa-ratio
    k_J_juv = 0.027, # juvenile and adult maturity maintenance rate constant; initially assumed equal to k_J_emb
    H_j1 = 1, # maturity at the start of metamorphosis (decline of feeding)
    H_p = 55., # maturity at puberty
    delta_E = 1., # energy density of E_mt, relative to remaining dry mass (not used)

    T_A = 8000., # Arrhenius temperature (K) 
    T_ref = 293.15, # reference 
    b_T = 40., # effect strength of temperature on resource allocation

    h_b_aq = 1e-3, 
    h_b_ter = 1e-3,

    watercontent_larvae = 0.9, 
    watercontent_juveniles = 0.75, 

    aux = ComponentVector(
        tau_R = 365., 
        h_S = 1e-3, 
        S_rel_crit = 1/3, 
        a_max = Truncated(Normal(3650, 365), 0, Inf)
    )
)

defaultparams() = ComponentVector(glb = glb, spc = spc)

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