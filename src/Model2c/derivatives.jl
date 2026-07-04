import ..AmphiDEB.Model2: food_dynamics_firstorder!
import ..AmphiDEB.Model2: y_T
import ..AmphiDEB.Model2: y_T_kap
import ..AmphiDEB.Model2: f_X
import ..AmphiDEB.Model2: calc_S_max
import ..AmphiDEB.Model2: embryo!
import ..AmphiDEB.Model2: larva!
#import ..AmphiDEB.Model2: metamorph!
#import ..AmphiDEB.Model2: juvenile!
import ..AmphiDEB.Model2: adult!
import ..AmphiDEB.Model2: sys_embryo!
import ..AmphiDEB.Model2: sys_larva!
#import ..AmphiDEB.Model2: sys_metamorph!
#import ..AmphiDEB.Model2: sys_juvenile!
import ..AmphiDEB.Model2: sys_adult!
#import ..AmphiDEB.Model2: sim_embryo
import ..AmphiDEB.Model2: sim_juvenile
import ..AmphiDEB.Model2: sim_adult
import ..AmphiDEB.Model2: isoutofdomain
import ..AmphiDEB.Model1: aquatic_tk_binary_loglogistic!, aquatic_td_binary_IA_loglogistic


function isoutofdomain(u, p, t)
    return u.ind.S < 0
end

"""
Smooth absolute function. 
Parameter λ controls the "smoothness".
"""
function sabs(x; λ = 1e-6)
    return sqrt(x^2 + λ)
end

"""
Smooth maximum function.
Parameter λ controls the "smoothness" via sabs.
"""
function smax(a, b; λ = 1e-6)
    return (a+b+sabs(a-b; λ = λ))/2
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
ODE component for the larval life stage.
"""
function larva!(du, u, p, t; 
    y_G = 1., y_GP = 1.,
    y_M = 1., y_MP = 1.,
    y_A = 1., y_AP = 1.
    )::Nothing

    # usimg max() in combination with isoutofdomain based on https://discourse.julialang.org/t/domainerror-while-solving-ode/53199/4, 
    u = max.(0, u)

    
    @unpack T_aq, V_patch_aq, food_dynamic = p.glb
    @unpack Z, dI_max_lrv, eta_IA, eta_AS_emb, eta_SA, kappa_emb, k_M_emb, k_M_Emt, k_J_emb, b_T, T_ref, T_A, K_X_lrv, gamma = p.ind
    @unpack X_aq = u.glb
    @unpack S, H, E_mt = u.ind
    
    aquatic_tk_binary_loglogistic!(du, u, p, t)

    y_G, y_M, y_A, _, _, y_κ_neg, y_κ_pos = aquatic_td_binary_IA_loglogistic(du, u, p, t)

    yT = y_T(T_A, T_ref, T_aq)
    fX = f_X(X_aq, V_patch_aq, K_X_lrv)

    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T_aq) * y_κ_neg * y_κ_pos

    S = max(0, S) # this needed for isoutofdomain() to work 

    dI = fX * dI_max_lrv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = (S * k_M_emb + E_mt * k_M_Emt) * y_M * y_MP * yT
    dJ = H * k_J_emb * y_M * y_MP * yT

    dS = Base.ifelse(
        kappa_T * dA >= dM, 
        y_G * y_GP * eta_AS_emb * (1 - gamma) * (kappa_T * dA - dM),
        -(dM / eta_SA - kappa_T * dA)
    )

    dE_mt = Base.ifelse(
        (kappa_T * dA) > dM, 
        eta_AS_emb * y_G * y_GP * gamma * (kappa_T * dA - dM),
        -(dM / eta_SA - (1 - gamma) * kappa_T * dA)
    )

    dE_mt_max = dE_mt
    dH = max(0, (1 - kappa_T) * dA - dJ)

    du.glb.X_aq -= (dI * p.glb.food_dynamic)
    du.glb.X_ter = 0.
    du.ind.X_emb = 0.
    du.ind.I = dI
    du.ind.A = dA
    du.ind.M = dM
    du.ind.J = dJ
    du.ind.S = dS
    du.ind.E_mt = dE_mt
    du.ind.E_mt_max = dE_mt_max
    du.ind.R = 0.
    du.ind.H = dH

    return nothing
end

"""
ODE component for the metamorph life stage.

Metamorphs consume the reserve buffer `E_mt` to pay costs for maintenance and maturation.
"""
function metamorph!(
    du, u, p, t;
    y_G = 1, y_GP = 1,
    y_M = 1, y_MP = 1,
    y_A = 1, y_AP = 1,
    )::Nothing

    # usimg max() in combination with isoutofdomain based on https://discourse.julialang.org/t/domainerror-while-solving-ode/53199/4, 
    u = max.(0, u)

    @unpack T_aq, V_patch_aq = p.glb
    @unpack dI_max_lrv, eta_IA, eta_AS_emb, eta_SA, kappa_emb, b_T, T_ref, K_X_lrv, k_M_emb, k_M_Emt, k_J_emb, k_J_juv, k_T, T_A = p.ind
    @unpack X_aq = u.glb
    @unpack E_mt, E_mt_max, S, H = u.ind

    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T_aq)
    yT = y_T(T_A, T_ref, T_aq)

    dI = 0.
    dA = 0.
    dM = ((S * k_M_emb) + (E_mt * k_M_Emt)) * y_M * y_MP * yT
    dJ = (H * k_J_emb * y_M * y_MP * yT) # maintenance costs for larval maturity

    dC = k_T * E_mt * yT # active mobilization flux
    dE_mt = -dC # TBD: add efficiency `/eta_EC` here? 

    dH = -dC  # larval maturity decreases according to a fixed rate
    dH = (1 - kappa_T) * dC - dJ # juvenile maturity is built during metamorphosis
    
    dS = Base.ifelse( 
        kappa_T * dC > dM, 
        (eta_AS_emb * y_G * y_GP * (kappa_T * dC - dM)), # structural growth and maintenance are fueled by E_mt according to κ-rule 
        -(dM / eta_SA - kappa_T * dC) # if the flux of E_mt is not sufficient, the shrinking equation applies
    )

    du.glb.X_aq -= (dI * p.glb.food_dynamic)
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
    du.ind.H = dH

    return nothing
end

"""
ODE component for the juvenile life stage.
"""
function juvenile!(
    du, u, p, t; 
    y_G = 1., y_GP = 1.,
    y_A = 1., y_AP = 1.,
    y_M = 1., y_MP = 1.,
    )::Nothing

    # usimg max() in combination with isoutofdomain based on https://discourse.julialang.org/t/domainerror-while-solving-ode/53199/4, 
    u = max.(0, u)

    @unpack T_ter, A_patch_ter, food_dynamic = p.glb
    @unpack X_ter = u.glb

    @unpack dI_max_juv, eta_IA, k_M_juv, k_M_Emt, k_J_emb, k_J_juv, eta_AS_juv, eta_SA, kappa_juv, b_T, T_ref, T_A, K_X_juv = p.ind
    @unpack S, H, E_mt = u.ind

    kappa_T = y_T_kap(kappa_juv, b_T, T_ref, T_ter)
    yT = y_T(T_A, T_ref, T_ter)
    fX = f_X(X_ter, A_patch_ter, K_X_juv) 

    dI = fX * dI_max_juv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = (S * k_M_juv + E_mt * k_M_Emt) * y_M * y_MP * yT
    dJ = H * k_J_emb * y_M * y_MP * yT 
    dS = Base.ifelse( 
        kappa_T * dA >= dM, 
        y_G * y_GP * eta_AS_juv * (kappa_T * dA - dM),
        -(dM / eta_SA - kappa_T * dA), 
    )    

    du.glb.X_aq = 0.
    du.glb.X_ter -= (food_dynamic * dI)
    du.ind.X_emb = 0.
    du.ind.E_mt = 0.
    du.ind.E_mt_max = 0.
    du.ind.I = dI 
    du.ind.A = dA
    du.ind.M = dM
    du.ind.J = dJ
    du.ind.S = dS
    du.ind.H = 0.
    du.ind.R = 0.

    return nothing
end


"""
ODE component for the adult life stage.
"""
function adult!(du, u, p, t;
    y_G = 1., y_GP = 1.,
    y_A = 1., y_AP = 1.,
    y_M = 1., y_MP = 1.,
    y_R = 1., y_MR = 1.
    )::Nothing

    # usimg max() in combination with isoutofdomain based on https://discourse.julialang.org/t/domainerror-while-solving-ode/53199/4, 
    u = max.(0, u)
    
    @unpack food_dynamic, T_ter, A_patch_ter = p.glb
    @unpack X_ter = u.glb

    @unpack dI_max_juv, eta_IA, k_M_juv, k_M_Emt, k_J_emb, k_J_juv, eta_AS_juv, eta_SA, eta_AR, kappa_juv, b_T, T_ref, T_A, K_X_juv = p.ind
    @unpack S, H, E_mt = u.ind

    kappa_T = y_T_kap(kappa_juv, b_T, T_ref, T_ter)
    yT = y_T(T_A, T_ref, T_ter)
    fX = f_X(X_ter, A_patch_ter, K_X_juv) 

    dI = fX * dI_max_juv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = (S * k_M_juv + E_mt * k_M_Emt) * y_M * y_MP * yT
    dJ = H * k_J_juv * y_M * y_MP * yT # adult pays maintenance for residual larval maturity
    dS = y_G * y_GP * eta_AS_juv * (kappa_T * dA - dM)

    dR =  eta_AR * y_R * y_MR * ((1 - kappa_T) * dA - dJ)
    
    du.glb.X_ter -= (dI * food_dynamic)
    du.ind.X_emb = 0.
    du.ind.E_mt = 0.
    du.ind.E_mt_max = 0.
    du.ind.I = dI 
    du.ind.A = dA
    du.ind.M = dM
    du.ind.J = dJ
    du.ind.S = dS
    du.ind.H = 0.
    du.ind.R = dR

    return nothing
end


function sys_metamorph!(du, u, p, t)
    food_dynamics_firstorder!(du, u, p, t)
    metamorph!(du, u, p, t)
    return nothing
end

"""
ODE system for juveniles including global component.
"""
function sys_juvenile!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du, u, p, t)
    juvenile!(du, u, p, t)
end

"""
ODE system for adults including global component.
"""
function sys_adult!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du, u, p, t)
    adult!(du, u, p, t)
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


"""
ODE system for all life stages, triggering parts of the ODE with conditional statements.
This is mostly useful for IBM simulations.
"""
function sys_individual_complete!(du, u, p, t)

    if u.ind.embryo>0 
        embryo!(du, u, p, t)
    end

    if u.ind.larva>0 
        larva!(du, u, p, t)
    end

    if u.ind.metamorph>0
        metamorph!(du, u, p, t)
    end

    if u.ind.juvenile>0
        juvenile!(du, u, p, t)
    end

    if u.ind.adult>0
        adult!(du, u, p, t)
    end
end