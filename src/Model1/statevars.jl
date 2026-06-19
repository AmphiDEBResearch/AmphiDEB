

const X_EMB_INT_REL = 1e-3

""""
    initialize_individual_statevars(p::ComponentVector; kwargs...)::ComponentVector

Initialize individual-level state variables for the AmhpiDEB model. 
Additional states can be added via kwargs. 
If existing states are provided via kwargs, these will be overwritten.
"""
function initialize_individual_statevars(p::ComponentVector; id = 1, cohort = 0)::ComponentVector

    return ComponentVector(
        embryo = 1.,
        larva = 0, # additional life stage: larva
        metamorph = 0, # additional life stage: metamorph
        juvenile = 0.,
        adult = 0.,

        X_emb = p.ind.X_emb_int, # initial mass of vitellus
        S = p.ind.X_emb_int * X_EMB_INT_REL, # initial structure is a small fraction of initial reserve // mass of vitellus
        H = 0., # maturity
        R = 0., # reproduction buffer
        f_X = 1., # scaled functional response 
        I_emb = 0., # ingestion from vitellus
        I = 0., # total ingestion
        A = 0., # assimilation
        M = 0., # somatic maintenance
        J = 0., # maturity maintenance 
        
        S_max = calc_S_max(p.ind.dI_max_emb, p.ind.eta_IA, p.ind.kappa_emb, p.ind.k_M_emb), # currently possible maximum stuructural mass

        E_mt = 1e-10, # metamorphic reserve
        E_mt_max = 1e-10, # maximum metamorphic reserve
        P_S = 0, # pathogen sporangia

        tktd = tktdstates_binary(),

        # auxiliary variables, only needed for population modelling
        aux = ComponentVector(
            id = id, 
            cohort = cohort,
            S_max_hist = p.ind.X_emb_int * X_EMB_INT_REL, # initial reference structure
            age = 0.,
            cause_of_death = 0.,
            time_since_last_repro = 0.,
            cum_repro = 0.,
            fX = 1.,
        )
    )
end

function tktdstates_binary()
    return ComponentVector(
        D_W1_G = 0.,
        D_W1_M = 0.,
        D_W1_A = 0.,
        D_W1_R = 0.,
        D_W1_H = 0.,
        D_W1_κ_pos = 0.,
        D_W1_κ_neg = 0.,

        D_W2_G = 0.,
        D_W2_M = 0.,
        D_W2_A = 0.,
        D_W2_R = 0.,
        D_W2_H = 0.,
        D_W2_κ_pos = 0.,
        D_W2_κ_neg = 0.,
    )
end

function initialize_global_statevars(p)
    return ComponentVector(
        X_aq = p.glb.dX_in_aq,
        X_ter = p.glb.dX_in_ter,
        C_W1 = 0., 
        C_W2 = 0.,
        P_Z = 0,
        N_emb = 0,
        N_lrv = 0,
        N_mt = 0,
        N_juv = 0,
        N_ad = 0,
        aging_mortality = 0, # cumulative mortalities
        starvation_mortality = 0,
        background_mortality = 0,
        S_mean = 0, # average structural mass
        S_sd = 0, # sd of of structural mass
    )
end


function initialize_statevars(p::ComponentVector)

    global_statevars = initialize_global_statevars(p)
    individual_statevars = initialize_individual_statevars(p)
    
    return ComponentVector(
        glb = global_statevars,
        ind = individual_statevars
    )
end