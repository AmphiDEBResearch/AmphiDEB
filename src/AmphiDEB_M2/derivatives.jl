# derivatives_M2.jl
# model version with alternative rules for climax: 
#   ingestion immediately goes to 0
#   we make the simplifying assumption that dS = 0 during metamorphosis 

"""
    AmphiDEB_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

ODE system with log-logistic toxicodynamics. 
"""
function M2_complete_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

    AmphiDEB_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    M2_individual_ODE_with_loglogistic_TD!(du, u, p, t)

    return nothing
end

function M2_individual_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

    TKTD_LL2!(du, u, p, t) # TKTD with mixtures assuming IA
    Pathogen_Infection!(du, u, p, t) # infection, release of zoospores and relative response to sporangia density
    M2_DEB!(du, u, p, t) # Amphibian DEB model

    return nothing
end


function M2_DEB!(du, u, p, t)::Nothing

    determine_life_stage!(du, u, p, t)
    Arrhenius!(du, u, p, t)
    eta_AS, kappa = life_stage_and_plasticity_effects!(du, u, p, t)

    M2_ingestion!(du, u, p, t)
    maintenance!(du, u, p, t)
    growth!(du, u, p, t, eta_AS, kappa)
    maturation!(du, u, p, t, kappa)
    M2_metamorphic_reserve!(du, u, p, t, eta_AS, kappa)

    reproduction!(du, u, p, t, kappa)

    return nothing
end


function M2_ingestion!(du, u, p, t)::Nothing

    u.ind.f_X = f_X(u.ind[:larva], u.ind[:metamorph], u.ind[:juvenile], u.ind[:adult], u.glb[:X], p.glb[:V_patch], p.ind[:K_X_lrv], p.ind[:K_X_juv])

    dI_emb = calc_dI_emb(u.ind[:embryo], u.ind[:S], p.ind[:dI_max_emb], u.ind[:y_T]) # ingestion by embryos
    dI_mt = 0. 
    dI_lrv = calc_dI_lrv(u.ind[:larva], u.ind[:f_X], p.ind[:dI_max_lrv], u.ind[:S], u.ind[:y_T]) # ingestion by larvae
    dI_juv_ad = calc_dI_juv(u.ind[:juvenile], u.ind[:adult], u.ind[:f_X], p.ind[:dI_max_juv], u.ind[:S], u.ind[:y_T])
    
    du.ind.I = dI(dI_emb, dI_mt, dI_lrv, dI_juv_ad)
    du.ind.X_emb = -dI_emb
    du.glb.X[1] -= (u.ind[:larva] + u.ind[:metamorph]) * du.ind.I
    du.glb.X[2] -= (u.ind[:juvenile] + u.ind[:adult]) * du.ind.I

    # assimilation flux
    du.ind.A = dA(du.ind[:I], p.ind[:eta_IA], u.ind[:y_j][3], u.ind[:y_jP][3])

    return nothing 
end

"""
    M2_calc_dE_mt_mt(
        dH::Float64,
        dJ::Float64,
        dM::Float64, 
        delta_E
        )::Float64

Calculation of reserve depletion (`dE_mt`) for metamorphs in model variant M2.

## Arguments

- dH: Maturation rate
- dJ: Maturity maintenance rate
- dM: Somatic maintenance rate
- delta_E: Energetic density of reserve, relative to remaining dry mass
"""
@inline function M2_calc_dE_mt_mt(
    dH::Float64,
    dJ::Float64,
    dM::Float64,
    delta_E::Float64
    )::Float64

    return -(dH + dJ + dM)/delta_E

end

function M2_metamorphic_reserve!(du, u, p, t, eta_AS::Float64, kappa::Float64)::Nothing
    
    dE_mt_lrv = calc_dE_mt_lrv(p.ind[:eta_AS_emb], u.ind[:y_j][1], u.ind[:y_jP][1], p.ind[:gamma], kappa, du.ind[:A], du.ind[:M], p.ind[:eta_SA], p.ind[:delta_E]) 
    dE_mt_mt = M2_calc_dE_mt_mt(du.ind[:H], du.ind[:J], du.ind[:M], p.ind[:delta_E])
    du.ind.E_mt = dE_mt(u.ind[:larva], dE_mt_lrv, u.ind[:metamorph], dE_mt_mt)
    du.ind.E_mt_max = dE_mt_max(u.ind[:larva], du.ind[:E_mt])

    return nothing
end