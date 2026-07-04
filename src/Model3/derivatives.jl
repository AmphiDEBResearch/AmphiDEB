"""
ODE component for food dynamics with simple first-order kinetics.
Disengaged for `p.glb.food_dynamic = 0.`
"""
function food_dynamics_firstorder!(du, u, p, t)::Nothing
    
    @unpack dX_in_aq, k_V_aq, dX_in_ter, k_V_ter, food_dynamic = p.glb
    @unpack X_aq, X_ter = u.glb

    du.glb.X_aq = food_dynamic * (dX_in_aq - k_V_aq * X_aq)
    du.glb.X_ter = food_dynamic * (dX_in_ter - k_V_ter * X_ter )

    return nothing 
end

"""
One-parameter Arrhenius temperature correction.
"""
@inline function y_T(
    T_A::Real,
    T_ref::Real,
    T::Real
    )::Real
    return exp((T_A / T_ref) - (T_A / T))
end

"""
Specific temperature correction on kappa according to Romoli et al.
"""
function y_T_kap(kappa, b_T, T_ref, T)
    1/(1 + (((1-kappa)/kappa) * exp(-b_T * ((T_ref - T)/T_ref))))
end

"""
Scaled functional response `f_X` based on Holling Type II functional response.
"""
@inline function f_X(
    X::Real,
    V_patch::Real,
    K_X::Real
    )::Real
    return (X / V_patch) / ((X / V_patch) + K_X)
end

@inline function calc_S_max(
    dI_max::Real, 
    eta_IA::Real, 
    kappa::Real, 
    k_M::Real
    )::Real
    return ((dI_max*eta_IA*kappa)/k_M)^3
end

"""
ODE component for the embryonic life stage.
"""
function embryo!(
    du, u, p, t; 
    y_G = 1., y_GP = 1., 
    y_M = 1., y_MP = 1.,
    y_A = 1., y_AP = 1.
    )::Nothing

    # usimg max() in combination with isoutofdomain based on https://discourse.julialang.org/t/domainerror-while-solving-ode/53199/4, 
    u = max.(0, u)

    @unpack T_aq, food_dynamic = p.glb
    @unpack dI_max_emb, eta_IA, k_M_emb, T_A, T_ref, k_J_emb, eta_AS_emb, kappa_emb = p.ind
    @unpack S, E_mt = u.ind

    yT = y_T(T_A, T_ref, T_aq)

    dI = S^(2/3) * dI_max_emb * yT
    dA = eta_IA * y_A * y_AP * dI
    dM = S * k_M_emb * y_M * y_MP * yT
    dJ = E_mt * k_J_emb * y_M * y_MP * yT
    dS = y_G * y_GP * eta_AS_emb * (kappa_emb * dA - dM)
    dE_mt =  smax(0, (1 - kappa_emb) * dA - dJ)
    
    du.ind.X_emb = -dI
    du.ind.I = dI
    du.ind.M = dM
    du.ind.J = dJ
    du.ind.S = dS
    du.ind.H = 0.
    du.ind.R = 0.
    du.ind.E_mt = dE_mt
    du.ind.E_mt_max = 0.

    return nothing
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
    @unpack Z, dI_max_lrv, eta_IA, eta_AS_emb, eta_SA, kappa_emb, k_M_emb, k_J_emb, b_T, T_ref, T_A, K_X_lrv = p.ind
    @unpack X_aq = u.glb
    @unpack S, E_mt = u.ind
    
    aquatic_tk_binary_loglogistic!(du, u, p, t)

    y_G, y_M, y_A, _, _, y_κ_neg, y_κ_pos = aquatic_td_binary_IA_loglogistic(du, u, p, t)

    yT = y_T(T_A, T_ref, T_aq)
    fX = f_X(X_aq, V_patch_aq, K_X_lrv)

    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T_aq) * y_κ_neg * y_κ_pos

    S = max(0, S)

    dI = fX * dI_max_lrv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M_emb * y_M * y_MP * yT
    dJ = E_mt * k_J_emb * y_M * y_MP * yT

    dS = Base.ifelse(
        kappa_T * dA >= dM, 
        y_G * y_GP * eta_AS_emb * (kappa_T * dA - dM),
        -(dM / eta_SA - kappa_T * dA)
    )
    
    dE_mt = smax(0, (1 - kappa_T) * dA - dJ)

    du.glb.X_aq -= (dI * p.glb.food_dynamic)
    du.glb.X_ter = 0.
    du.ind.X_emb = 0.
    du.ind.I = dI
    du.ind.A = dA
    du.ind.M = dM
    du.ind.J = dJ
    du.ind.S = dS
    du.ind.E_mt = dE_mt
    du.ind.E_mt_max = dE_mt
    du.ind.R = 0.
    du.ind.H = 0.

    return nothing
end

"""
ODE component for the metamorph life stage.
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
    @unpack dI_max_lrv, eta_IA, eta_AS_emb, kappa_emb, b_T, T_ref, K_X_lrv, k_M_emb, k_J_emb, k_J_juv, k_C, T_A = p.ind
    @unpack X_aq = u.glb
    @unpack E_mt, E_mt_max, H,  S = u.ind

    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T_aq)
    yT = y_T(T_A, T_ref, T_aq)
    fX = f_X(X_aq, V_patch_aq, K_X_lrv)

    dI_max_t = dI_max_lrv * (E_mt / E_mt_max)
    dI = fX * dI_max_t * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M_emb * y_M * y_MP * yT
    dJ = ((E_mt * k_J_emb) + (H * k_J_juv)) * y_M * y_MP * yT
    dH = max(0, ((1 - kappa_T) * dM) / kappa_T - dJ)
    
    dE_mt = -k_C * E_mt # reserve buffer decreases
    dC = eta_IA * (dI + dE_mt) * yT # mobliziation flux equals ingestion + reserve buffer mobilization
    dS = eta_AS_emb * y_G * ((kappa_T * dC) - dM)
    dH = ((1-kappa_T) * dC) - dJ

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
    du.ind.E_mt_max = 0.

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

    @unpack dI_max_juv, eta_IA, k_M_juv, k_J_juv, eta_AS_juv, eta_SA, kappa_juv, b_T, T_ref, T_A, K_X_juv = p.ind
    @unpack S, H = u.ind

    kappa_T = y_T_kap(kappa_juv, b_T, T_ref, T_ter)
    yT = y_T(T_A, T_ref, T_ter)
    fX = f_X(X_ter, A_patch_ter, K_X_juv) 

    dI = fX * dI_max_juv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M_juv * y_M * y_MP * yT
    dJ = H * k_J_juv * y_M * y_MP * yT 
    dS = Base.ifelse( 
        kappa_T * dA >= dM, 
        y_G * y_GP * eta_AS_juv * (kappa_T * dA - dM),
        -(dM / eta_SA - kappa_T * dA), 
    )    
    dH =  max(0, (1 - kappa_T) * dA - dJ)

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
    du.ind.H = dH
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

    @unpack dI_max_juv, eta_IA, k_M_juv, k_J_juv, eta_AS_juv, eta_SA, eta_AR, kappa_juv, b_T, T_ref, T_A, K_X_juv = p.ind
    @unpack S, H = u.ind

    kappa_T = y_T_kap(kappa_juv, b_T, T_ref, T_ter)
    yT = y_T(T_A, T_ref, T_ter)
    fX = f_X(X_ter, A_patch_ter, K_X_juv) 

    dI = fX * dI_max_juv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M_juv * y_M * y_MP * yT
    dJ = H * k_J_juv * y_M * y_MP * yT 
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

minimal_TK(k_D, C_W, D_W) = k_D * (C_W - D_W)

function aquatic_tk_binary_loglogistic!(du, u, p, t)::Nothing

    p_tktd = p.ind.tktd
    u_tk = u.ind.tktd

    du.ind.tktd.D_W1_G = minimal_TK(p_tktd.k_D1_G, p.glb.C_W1, u_tk.D_W1_G)
    du.ind.tktd.D_W1_M = minimal_TK(p_tktd.k_D1_M, p.glb.C_W1, u_tk.D_W1_M)
    du.ind.tktd.D_W1_A = minimal_TK(p_tktd.k_D1_A, p.glb.C_W1, u_tk.D_W1_A)
    du.ind.tktd.D_W1_R = minimal_TK(p_tktd.k_D1_R, p.glb.C_W1, u_tk.D_W1_R)
    du.ind.tktd.D_W1_H = minimal_TK(p_tktd.k_D1_H, p.glb.C_W1, u_tk.D_W1_H)
    du.ind.tktd.D_W1_κ_pos = minimal_TK(p_tktd.k_D1_κ_pos, p.glb.C_W1, u_tk.D_W1_κ_neg)
    du.ind.tktd.D_W1_κ_neg = minimal_TK(p_tktd.k_D1_κ_neg, p.glb.C_W1, u_tk.D_W1_κ_pos)
    
    du.ind.tktd.D_W2_G = minimal_TK(p_tktd.k_D2_G, p.glb.C_W2, u_tk.D_W2_G)
    du.ind.tktd.D_W2_M = minimal_TK(p_tktd.k_D2_M, p.glb.C_W2, u_tk.D_W2_M)
    du.ind.tktd.D_W2_A = minimal_TK(p_tktd.k_D2_A, p.glb.C_W2, u_tk.D_W2_A)
    du.ind.tktd.D_W2_R = minimal_TK(p_tktd.k_D2_R, p.glb.C_W2, u_tk.D_W2_R)
    du.ind.tktd.D_W2_H = minimal_TK(p_tktd.k_D2_H, p.glb.C_W2, u_tk.D_W2_H)
    du.ind.tktd.D_W2_κ_pos = minimal_TK(p_tktd.k_D2_κ_pos, p.glb.C_W2, u_tk.D_W2_κ_neg)
    du.ind.tktd.D_W2_κ_neg = minimal_TK(p_tktd.k_D2_κ_neg, p.glb.C_W2, u_tk.D_W2_κ_pos)
    
    return nothing
end

LL2neg(D, e, b) = 1/(1 + (D/e)^b)
LL2pos(D, e, b) = 1 - log(LL2neg(D, e, b))

function aquatic_td_binary_IA_loglogistic(du, u, p, t)

    u_tk = u.ind.tktd 
    p_tktd = p.ind.tktd

    y_G1 = LL2neg(u_tk.D_W1_G, p_tktd.e1_G, p_tktd.b1_G)
    y_G2 = LL2neg(u_tk.D_W2_G, p_tktd.e2_G, p_tktd.b2_G)
    y_G = y_G1 * y_G2

    y_M1 = LL2pos(u_tk.D_W1_M, p_tktd.e1_M, p_tktd.b1_M)
    y_M2 = LL2pos(u_tk.D_W2_M, p_tktd.e2_M, p_tktd.b2_M)
    y_M = y_M1 * y_M2

    y_A1 = LL2neg(u_tk.D_W1_A, p_tktd.e1_A, p_tktd.b1_A)
    y_A2 = LL2neg(u_tk.D_W2_A, p_tktd.e2_A, p_tktd.b2_A)
    y_A = y_A1 * y_A2

    y_R1 = LL2neg(u_tk.D_W1_R, p_tktd.e1_R, p_tktd.b1_R)
    y_R2 = LL2neg(u_tk.D_W2_R, p_tktd.e2_R, p_tktd.b2_R)
    y_R = y_R1 * y_R2

    y_H1 = LL2pos(u_tk.D_W1_H, p_tktd.e1_H, p_tktd.b1_H)
    y_H2 = LL2pos(u_tk.D_W2_H, p_tktd.e2_H, p_tktd.b2_H)
    y_H = y_H1 * y_H2

    y_κ_pos_1 = LL2neg(u_tk.D_W1_κ_pos, p_tktd.e1_κ_pos, p_tktd.b1_κ_pos)
    y_κ_pos_2 = LL2neg(u_tk.D_W2_κ_pos, p_tktd.e2_κ_pos, p_tktd.b2_κ_pos)
    y_κ_pos = y_κ_pos_1 * y_κ_pos_2

    y_κ_neg_1 = LL2neg(u_tk.D_W1_κ_neg, p_tktd.e1_κ_neg, p_tktd.b1_κ_neg)
    y_κ_neg_2 = LL2neg(u_tk.D_W2_κ_neg, p_tktd.e2_κ_neg, p_tktd.b2_κ_neg)
    y_κ_neg = y_κ_neg_1 * y_κ_neg_2

    return (y_G, y_M, y_A, y_R, y_H, y_κ_neg, y_κ_pos)
end

"""
ODE system for embryos including global component.
"""
function sys_embryo!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du, u, p, t)
    embryo!(du, u, p, t)
end

"""
ODE system for larvae including global component.
"""
function sys_larva!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du, u, p, t)
    larva!(du, u, p, t)
end

"""
ODE system for metamorphs including global component.
"""
function sys_metamorph!(du, u, p, t)::Nothing
    food_dynamics_firstorder!(du, u, p, t)
    metamorph!(du, u, p, t)
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
Complete ODE with if/else for triggering life stage-specific derivatives. 
Suitable for use in IBM context.
"""
function individual_ODE!(du, u, p, t)::Nothing

    if u.ind.embryo > 0 
        embryo!(du, u, p, t)
    end

    if u.ind.larva > 0
        larva!(du, u, p, t)
    end

    if u.ind.metamorph > 0
        metamorph!(du, u, p, t)
    end

    if u.ind.juvenile > 0
        juvenile!(du, u, p, t)
    end

    if u.ind.adult > 0 
        adult!(du, u, p, t)
    end

    return nothing
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
Simulate model variant without κ and first-order mobilization of E_mt.
This is inspired by the DEBlipis model (Martin et al. 2017), in the sense that the maturity takes over the role of a buffer compartment (E_mt).
During metamorphosis, we assume that E_mt is mobilized according to first-order kinetics and that size-specific ingestion rate dI_max drop at the same relative rate,
leading to exponential decline of both E_mt and dI_max. 

The combined mobilized flux E_mt + dI_max is disributed according to the κ-rule.
Maturity is built up during metamorphosis according to standard DEB(kiss) dynamics.

## References

Martin, B. T., Heintz, R., Danner, E. M., & Nisbet, R. M. (2017). Integrating lipid storage into general representations of fish energetics. Journal of Animal Ecology, 86(4), 812-825.
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

