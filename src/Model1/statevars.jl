

const X_EMB_INT_REL = 1e-3

""""
    initialize_individual_statevars(p::ComponentVector; kwargs...)::ComponentVector

Initialize individual-level state variables for the AmhpiDEB model. 
Additional states can be added via kwargs. 
If existing states are provided via kwargs, these will be overwritten.

**IMPORTANT NOTE FOR SIMULATING MIXTURES**: 
You can simulate an arbitrary number of chemical stressors, 
but it is currently not possible to dynamically change the shape of vectors and matrices contained in a `ComponentVector`. 
In practice, this means: If you want to simulate mixtures, you cannot simply provide the parameters as kwargs to this function, 
but the entire `ComponentVector` has to be re-defined. 
You can do so by copy-pasting the definition body of this function and changing the shape of `y_j`. 
The same is true for the parameter vector, where the shape of all TKTD parameters has to be adjusted. 
An example is given in the unit tests of the `EcotoxSystems` package: https://github.com/SimonHansul/EcotoxSystems.jl/blob/main/test/test05_mixtures.jl
For application to the AmphiDEB model, we have to take into account that it has an additional PMoA, and therefore an additional column in the sublethal TKTD parameters.
"""
function initialize_individual_statevars(p::ComponentVector)::ComponentVector

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

        # auxiliary variables, only needed for population modelling
        aux = ComponentVector(
            S_max_hist = p.ind.X_emb_int * X_EMB_INT_REL, # initial reference structure
            age = 0.,
            cause_of_death = 0.,
            time_since_last_repro = 0.,
            cum_repro = 0.,
            fX = 1.,
        )
    )
end


function initialize_global_statevars(p)
    return ComponentVector(
        N = 0, 
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
        aging_mortality = 0, # cumulative mortalities (only tracked in population modelling)
        starvation_mortality = 0,
        GUTS_mortality = 0,
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