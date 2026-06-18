import ..AmphiDEB.Model1: food_dynamics_firstorder!
import ..AmphiDEB.Model1: y_T
import ..AmphiDEB.Model1: y_T_kap
import ..AmphiDEB.Model1: f_X
import ..AmphiDEB.Model1: calc_S_max
import ..AmphiDEB.Model1: embryo!
import ..AmphiDEB.Model1: larva!
import ..AmphiDEB.Model1: juvenile!
import ..AmphiDEB.Model1: adult!

"""
ODE component for the metamorph life stage. 

Model 2: We enforce `dI = 0` and assume that `E_mt` is mobilized with a constant first-order rate `k_C`. Mobilized `E_mt` is distributed according to the κ-rule.
"""
function metamorph!(
    du, u, p, t;
    y_G = 1, y_GP = 1,
    y_M = 1, y_MP = 1,
    y_A = 1, y_AP = 1,
    )::Nothing

    # usimg max() in combination with isoutofdomain based on https://discourse.julialang.org/t/domainerror-while-solving-ode/53199/4, 
    u = max.(0, u)

    @unpack T_aq, V_patch_aq, food_dynamic = p.glb
    @unpack k_C, eta_IA, eta_AS_emb, kappa_emb, b_T, T_ref, K_X_lrv, k_M_emb, k_J_emb, delta_E, delta_k_M_mt, T_A = p.ind
    @unpack X_aq = u.glb
    @unpack E_mt, E_mt_max, S, H = u.ind

    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T_aq) # temperature-dependent κ
    k_M = k_M_emb * delta_k_M_mt 
    yT = y_T(T_A, T_ref, T_aq)

    pE_mt = k_C * E_mt

    dI = 0. 
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M * y_M * y_MP * yT
    dJ = H * k_J_emb * y_M * y_MP * yT

    dS = eta_AS_emb * y_G * y_GP * (kappa_T * pE_mt - dM)
    dH = (1 - kappa_T) * pE_mt - dJ
    
    dE_mt = -pE_mt
    
    du.glb.X_aq -= (food_dynamic * dI)
    du.glb.X_ter = 0.
    du.ind.X_emb = 0.
    du.ind.I = dI
    du.ind.A = dA
    du.ind.M = dM
    du.ind.S = dS
    du.ind.H = dH
    du.ind.J = dJ
    du.ind.R = 0.
    du.ind.E_mt = dE_mt
    du.ind.E_mt_max = 0.

    return nothing
end

"""
ODE system for embryos including global component.
"""
function sys_embryo!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du.glb, u.glb, p.glb, t)
    embryo!(du, u, p, t)
end

"""
ODE system for larvae including global component.
"""
function sys_larva!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du.glb, u.glb, p.glb, t)
    larva!(du, u, p, t)
end

"""
ODE system for metamorphs including global component.
"""
function sys_metamorph!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du.glb, u.glb, p.glb, t)
    metamorph!(du, u, p, t)
end

"""
ODE system for juveniles including global component.
"""
function sys_juvenile!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du.glb, u.glb, p.glb, t)
    juvenile!(du, u, p, t)
end

"""
ODE system for adults including global component.
"""
function sys_adult!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du.glb, u.glb, p.glb, t)
    adult!(du, u, p, t)
end


function isoutofdomain(u, p, t)
    return u.ind.S < 0
end

"""
Simulate embryo from initialization to birth.
""" 
function sim_embryo(p; saveat = [], alg = Rodas5P(), return_sol = false)

    p_ind = generate_individual_params(p)
    u0 = initialize_statevars(p_ind)
    tspan = (0,p.glb.t_max)

    prob = ODEProblem(sys_embryo!, u0, tspan, p_ind)
    sol = solve(prob, callback = birth_terminal, saveat = saveat, alg = alg, isoutofdomain = isoutofdomain)

    if return_sol
        return sol, sol.u[end], p_ind
    end

    return EcotoxSystems.sol_to_df(sol), sol.u[end], p_ind
end

"""
Simulate larva from birth to metamorphosis.
"""
function sim_larva(p_ind, u0; saveat = [], alg = Rodas5P())

    tspan = (0,p_ind.glb.t_max)
    
    prob = ODEProblem(sys_larva!, u0, tspan, p_ind)
    sol = solve(prob, callback = metamorphosis_terminal, saveat = saveat, alg = alg, isoutofdomain = isoutofdomain)
    
    return EcotoxSystems.sol_to_df(sol), sol.u[end], p_ind
end

"""
Simulate metamorph from metamorphosis (Gosner 42) to froglet emergence (Gosner 46)-
"""
function sim_metamorph(p_ind, u0; saveat = [], alg = Rodas5P())

    tspan = (0,p_ind.glb.t_max)

    prob = ODEProblem(sys_metamorph!, u0, tspan, p_ind)
    sol = solve(prob, callback = froglet_emergence_terminal, saveat = saveat, alg = alg, isoutofdomain = isoutofdomain)

    return EcotoxSystems.sol_to_df(sol), sol.u[end], p_ind 
end

"""
Simulate juvenile from froglet emergence (Gosner 46) to puberty.
"""
function sim_juvenile(p_ind, u0; saveat = [], alg = Rodas5P())

    tspan = (0,p_ind.glb.t_max)

    prob = ODEProblem(sys_juvenile!, u0, tspan, p_ind, saveat = saveat, alg = alg, isoutofdomain = isoutofdomain)
    sol = solve(prob, callback = puberty_terminal)

    return EcotoxSystems.sol_to_df(sol), sol.u[end], p_ind
end

"""
Simulate adult from puberty to pre-defined maximum time `p.glb.t_max`.
"""
function sim_adult(p_ind, u0; saveat = [], alg = Rodas5P())

    tspan = (0,p_ind.glb.t_max)

    prob = ODEProblem(sys_adult!, u0, tspan, p_ind)
    sol = solve(prob, alg = alg, saveat = saveat, isoutofdomain = isoutofdomain)

    return EcotoxSystems.sol_to_df(sol), sol.u[end], p_ind
end


"""
Simulate all life stages consecutively as separate ODE systems.
"""
function sim_all(p; kwargs...)

    sim_emb, u0lrv, p_ind = sim_embryo(p; kwargs...)

    if sim_emb.t[end] >= p.glb.t_max
        return sim_emb
    end
    
    sim_lrv, u0mt, p_ind = sim_larva(p_ind, u0lrv; kwargs...)
    sim_lrv[!,:t] = sim_lrv.t .+ sim_emb.t[end]

    if sim_lrv.t[end] >= p.glb.t_max
        return vcat(sim_emb, sim_lrv)
    end

    sim_mt, u0juv, p_ind = sim_metamorph(p_ind, u0mt; kwargs...)
    sim_mt[!,:t] = sim_mt.t .+ sim_lrv.t[end]

    if sim_mt.t[end] >= p.glb.t_max
        return vcat(sim_emb, sim_lrv, sim_mt)
    end

    sim_juv, u0ad, p_ind = sim_juvenile(p_ind, u0juv; kwargs...)
    sim_juv[!,:t] = sim_juv.t .+ sim_mt.t[end]
    
    if sim_juv.t[end] >= p.glb.t_max
        return vcat(sim_emb, sim_lrv, sim_mt, sim_juv)
    end

    sim_ad, uend, p_ind = sim_adult(p_ind, u0ad; kwargs...)
    sim_ad[!,:t] = sim_ad.t .+ sim_juv.t[end]
    
    return vcat(sim_emb, sim_lrv, sim_mt, sim_juv, sim_ad)
end


#function sys_complete!(du, u, p, t)
#    food_dynamics_firstorder!(du.glb, u.glb, p.glb, t)
#
#    if u.ind.is_embryo>0 
#        embryo!(du, u, p, t)
#    end
#
#    if u.ind.is_larva>0 
#        larva!(du, u, p, t)
#    end
#
#    if u.ind.is_metamorph>0
#        metamorph!(du, u, p, t)
#    end
#end