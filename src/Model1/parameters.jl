
"""
Global parameters with defaults.
"""
glb = ComponentVector(
    t_max = 56., # 8 - week simulation
    N0 = 1., # start with single value [] - only possible setting for ODE_simulator
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
    pathogen_inoculation_time = 30., # time-point of pathogen inoculation [t]
    medium_renewals = [0.] # time-points on which media renewals occur; results in removal of spores
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
    delta_k_M_mt = 1., # metamorph somatic maintenance rate constant, relative to that of embryos andd larvae
    k_J_emb = 0.027, # embryonic to metamorph maturity maintenance rate constant; initially calculated based on (1-kappa)/kappa-ratio
    k_J_juv = 0.027, # juvenile and adult maturity maintenance rate constant; initially assumed equal to k_J_emb
    H_j1 = 1, # maturity at the start of metamorphosis (decline of feeding)
    H_p = 55., # maturity at puberty
    delta_E = 1., # energy density of E_mt, relative to remaining dry mass

    T_A = 8000., # Arrhenius temperature (K) 
    T_ref = 293.15, # reference 
    b_T = 40., # effect strength of temperature on resource allocation

    #=
    TKTD parameters    
    =#

    # TK feedbacks (1.0/0.0 for on/off) 
    # only including growth feedback for now
    fb_G = 0., 

    h_b = 0., # background mortality
   
    # toxicity parameters substance 1

    k_D1_G = 0., # parameters for sublethal effects
    k_D1_M = 0., 
    k_D1_A = 0.,
    k_D1_R = 0.,
    k_D1_Hneg = 0.,
    k_D1_Hpos = 0.,
    k_D1_KAPneg = 0.,

    b1_G = 2., 
    b1_M = 2., 
    b1_A = 2.,
    b1_R = 2.,
    b1_Hneg = 2., # decrease in maturity threshold for metamorphosis
    b1_Hpos = 2., # increase in maturity threshold for metamorphosis
    b1_KAPneg = 2., # decrease in kappa

    e1_G = 1.,
    e1_M = 1., 
    e1_A = 1.,
    e1_R = 1.,
    e1_Hneg = 1.,
    e1_Hpos = 1.,
    e1_KAPneg = 1.,

    k_D1_h = 0., 

    E_h1 = 1., # parameters for lethal effecs

    b_h1 = 2., 
    C_h1 = 1., # proportionality constant to convert relative response to hazard rate - may optionally be used to link sublethal effects during larval period to lethal effects after metamorphosis

    
    # toxicity parameters substance 2

    k_D2_G = 0., # parameters for sublethal effects
    k_D2_M = 0., 
    k_D2_A = 0.,
    k_D2_R = 0.,
    k_D2_Hneg = 0.,
    k_D2_Hpos = 0.,
    k_D2_KAPneg = 0.,

    b2_G = 2., 
    b2_M = 2., 
    b2_A = 2.,
    b2_R = 2.,
    b2_Hneg = 2., # decrease in maturity threshold for metamorphosis
    b2_Hpos = 2., # increase in maturity threshold for metamorphosis
    b2_KAPneg = 2., # decrease in kappa

    e2_G = 1.,
    e2_M = 1., 
    e2_A = 1.,
    e2_R = 1.,
    e2_Hneg = 1.,
    e2_Hpos = 1.,
    e2_KAPneg = 1.,

    k_D2_h = 0.,  # parameters for lethal effecs 
    E_h2 = 1.,
    b_h2 = 2., 
    C_h2 = 1., # proportionality constant to convert relative response to hazard rate - may optionally be used to link sublethal effects during larval period to lethal effects after metamorphosis


    # following parameters are curently only used in an individual-based context, but could find application in the pure-ODE implementation 
    # for example by triggering emptying of the reproduction buffer through callbacks

    S_rel_crit = 0.66, # initial guess on how much body mass can be lost - X. laevis can lose up to 45% of mass, but accompanied with considerable starvation mortality (Merkle & Hanke (1987) Comp. Biochem. Physiol.)
    h_S = 0.6, # hazard rate below critical mass - ca. 50% daily survival probability ()
    a_max = truncated(Normal(15 * 365, 1.5 * 365), 0, Inf), # maximum age [d]
    tau_R = 365., # reproduction period [d]
    
    #=
    Pathogen dynamics and effect parameters
    =#

    Chi = LogNormal(log(1)+1^2, 1), # killing rate modifier, log-normal distribution with mode 1 and sigma 1, (mu = ln(mode)+sigma^2).
    e_PG = 1.,
    e_PM = 1.,
    e_PA = 1.,
    e_PR = 1.,
    b_PG = 2.,
    b_PM = 2.,
    b_PA = 2.,
    b_PR = 2.,
)

params = ComponentVector(glb = glb, spc = spc)

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