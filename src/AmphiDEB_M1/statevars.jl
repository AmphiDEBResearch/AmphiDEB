

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
function initialize_individual_statevars(
    p::ComponentVector; 
    id = 1., 
    cohort = 0.,
    kwargs...)::ComponentVector

    return ComponentVector(
        #### state variables of the default DEB model #### 

        embryo = 1.,
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
        
        S_max = calc_S_max(p.ind[:dI_max_emb], p.ind[:eta_IA], p.ind[:kappa_emb], p.ind[:k_M_emb]), # currently possible maximum stuructural mass
        chi_G = 1., # TK feedback growth

        D_j = EcotoxSystems.constrmmat(p.ind.KD), # sublethal damage per stressor and PMoA; constrmmat constructs a MMatrix from the shape of the input
        D_h = EcotoxSystems.constrmvec(p.ind.KD_h), # lethal damage per stressor; constrmvec constructs a MVector from the shape of the 

        y_T = 1.,
        
        y_j = @MMatrix([1. 1. 1. 1. 1. 1. 1.]),  # relative response per stressor and pmoa
        h_z = 0., # hazard rate caused by chemical stressors
        S_z = 1., # chemical-related survival probability

        # these are curently only needed in the IBM version, 
        # but may find application in the pure-ODE implementation 

        S_max_hist = p.ind.X_emb_int * X_EMB_INT_REL, # initial reference structure
        id = id, 
        cohort = cohort,
        age = 0.,
        cause_of_death = 0.,
        time_since_last_repro = 0.,
        cum_repro = 0.,
        
        #### additional state variables ####

        larva = 0, # additional life stage: larva
        metamorph = 0, # additional life stage: metamorph
        E_mt = 1e-10, # metamorphic reserve
        E_mt_max = 1e-10, # maximum metamorphic reserve
        P_S = 0, # pathogen sporangia
        # this should use EcotoxSystems.constrmvec to construct a mutable static array
        # it doesn't work right now, not sure why, so I use regular arrays 
        # maybe slows down simulation a little compared to static array, but otherwise the same
        y_jP = @MMatrix([1. 1. 1. 1.]), #EcotoxSystems.constrmvec(Vector(p.ind.e_P), fillval = 1) # response to pathogen per PMoA
        kwargs...
    )
end


function initialize_global_statevars(p)
    return ComponentVector(
        EcotoxSystems.initialize_global_statevars(p);
        aging_mortality = 0,
        starvation_mortality = 0,
        GUTS_mortality = 0,
        P_Z = 0,
        N_emb = 0,
        N_lrv = 0,
        N_mt = 0,
        N_juv = 0,
        N_ad = 0
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