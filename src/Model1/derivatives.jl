birth_condition(u, t, integrator) = u.ind.X_emb
function birth_affect!(integrator)
    integrator.u.ind.is_embryo = 0.
    integrator.u.ind.is_larva = 1.
    integrator.u.ind.is_metamorph = 0.
    integrator.u.ind.is_juvenile = 0.
    integrator.u.ind.is_adult = 0.
end

birth = ContinuousCallback(
    birth_condition, 
    birth_affect!, 
    nothing # we need `neg_affect! = nothing` to tell the solver that the affect should only occur for upcrossings (condition function switches from negative to positive) 
    )

function metamorphosis_condition(u, t, integrator)
    u.ind.H - integrator.p.ind.H_j1 # NOTE: for chemical effects, make sure that the effect is applied on H_j1 here
end

function metamorphosis_affect!(integrator)
    integrator.u.ind.is_embryo = 0.
    integrator.u.ind.is_larva = 0.
    integrator.u.ind.is_metamorph = 1.
    integrator.u.ind.is_juvenile = 0.
    integrator.u.ind.is_adult = 0.
end

metamorphosis = ContinuousCallback(
    metamorphosis_condition, 
    metamorphosis_affect!
)

function froglet_emergence_condition(u, t, integrator)
    (u.ind.E_mt <= 0) && (u.ind.H >= integrator.p.ind.H_j1) && (u.ind.H < integrator.p.ind.H_p)
end

function froglet_emergence_affect!(integrator)
    integrator.u.ind.is_embryo = 0.
    integrator.u.ind.is_larva = 0.
    integrator.u.ind.is_metamorph = 0.
    integrator.u.ind.is_juvenile = 1.
    integrator.u.ind.is_adult = 0.
end

froglet_emergence = DiscreteCallback(
    froglet_emergence_condition, 
    froglet_emergence_affect!
)

puberty_condition(u, t, integrator) = u.ind.H - integrator.p.ind.H_p
function puberty_affect!(integrator)
    integrator.u.ind.is_embryo = 0.
    integrator.u.ind.is_larva = 0.
    integrator.u.ind.is_metamorph = 0.
    integrator.u.ind.is_juvenile = 0.
    integrator.u.ind.is_adult = 1.
end
puberty = ContinuousCallback(puberty_condition, puberty_affect!, nothing)

callback_set = CallbackSet(birth, metamorphosis, froglet_emergence, puberty)

function global!(du, u, p, t)::Nothing
    
    @unpack dX_in_aq, k_V_aq, dX_in_ter, k_V_ter = p
    @unpack X_aq, X_ter = u

    du.X_aq = dX_in_aq - k_V_aq * X_aq 
    du.X_ter = dX_in_ter - k_V_ter * X_ter 

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

function embryo!(
    du, u, p, t; 
    y_G = 1., y_GP = 1., 
    y_M = 1., y_MP = 1.,
    y_A = 1., y_AP = 1.
    )::Nothing

    @unpack dI_max_emb, eta_IA, k_M_emb, T_A, T_ref, k_J_emb, eta_AS_emb, kappa_emb = p.ind
    @unpack T_aq = p.glb
    @unpack S, H = u.ind

    yT = y_T(T_A, T_ref, T_aq)

    dI = S^(2/3) * dI_max_emb * yT
    dA = eta_IA * y_A * y_AP * dI
    dM = S * k_M_emb * y_M * y_MP * yT
    dJ = H * k_J_emb * y_M * y_MP * yT
    dS = y_G * y_GP * eta_AS_emb * (kappa_emb * dA - dM)
    dH =  max(0, (1 - kappa_emb) * dA - dJ)
    
    du.ind.X_emb = -dI
    du.ind.I = dI
    du.ind.M = dM
    du.ind.J = dJ
    du.ind.S = dS
    du.ind.R = 0.
    du.ind.H = dH
    du.ind.E_mt = 0.
    du.ind.E_mt_max = 0.

    return nothing
end


function larva!(du, u, p, t; 
    y_G = 1., y_GP = 1.,
    y_M = 1., y_MP = 1.,
    y_A = 1., y_AP = 1.
    )::Nothing

    @unpack dI_max_lrv, eta_IA, eta_AS_emb, eta_SA, kappa_emb, k_M_emb, k_J_emb, b_T, T_ref, T_A, K_X_lrv, gamma, delta_E = p.ind
    @unpack T_aq, V_patch_aq = p.glb
    @unpack X_aq = u.glb
    @unpack S, H = u.ind

    yT = y_T(T_A, T_ref, T_aq)
    fX = f_X(X_aq, V_patch_aq, K_X_lrv)
    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T_aq)

    dI = fX * dI_max_lrv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M_emb * y_M * y_MP * yT
    dJ = H * k_J_emb * y_M * y_MP * yT

    dS = Base.ifelse(
        kappa_T * dA >= dM, 
        y_G * y_GP * eta_AS_emb * (kappa_T * dA - dM),
        -(dM / eta_SA - kappa_T * dA)
    )

    dE_mt = Base.ifelse(
        (kappa_T * dA) > dM, 
        eta_AS_emb * y_G * y_GP * gamma * (kappa_T * dA - dM),
        -(dM / eta_SA - (1 - gamma) * kappa_T * dA)/delta_E
    )

    dE_mt_max = dE_mt
    dH = max(0, (1 - kappa_T) * dA - dJ)

    du.glb.X_aq = -dI
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


function metamorph!(
    du, u, p, t;
    y_G = 1, y_GP = 1,
    y_M = 1, y_MP = 1,
    y_A = 1, y_AP = 1,
    )::Nothing

    @unpack T_aq, V_patch_aq = p.glb
    @unpack dI_max_lrv, eta_IA, eta_AS_emb, kappa_emb, b_T, T_ref, K_X_lrv, k_M_emb, k_J_emb, delta_E, delta_k_M_mt, T_A = p.ind
    @unpack X_aq = u.glb
    @unpack E_mt, E_mt_max, S, H = u.ind

    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T_aq)
    k_M = k_M_emb * delta_k_M_mt 
    yT = y_T(T_A, T_ref, T_aq)
    fX = f_X(X_aq, V_patch_aq, K_X_lrv)

    S = max(1e-10, S)

    dI = fX * dI_max_lrv * (E_mt/E_mt_max) * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M * y_M * y_MP * yT
    dJ = H * k_J_emb * y_M * y_MP * yT
    dS = eta_AS_emb * y_G * y_GP * dA
    dH = max(0, ((1 - kappa_T) * dM) / kappa_T - dJ)
    dE_mt = -(dH + dJ + dM)/delta_E
    
    @show fX S dS eta_AS_emb dA

    du.glb.X_aq = -dI
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


function juvenile!(du, u, p, t)::Nothing

    return nothing
end

function adult!(du, u, p, t)::Nothing

    return nothing
end

function sim_embryo(p)

    p_ind = Model1.generate_individual_params(p)
    u0 = Model1.initialize_statevars(p_ind)
    tspan = (0,p.glb.t_max)

    prob = ODEProblem(embryo!, u0, tspan, p_ind)
    sol = solve(prob, callback = Model1.callback_set)

    return EcotoxSystems.sol_to_df(sol), sol.u[end], p_ind
end

function sim_larva(p_ind, u0)


end

function du!(du, u, p, t)
    global!(du.glb, u.glb, p.glb, t)

    if u.ind.is_embryo>0 
        embryo!(du, u, p, t)
    end

    if u.ind.is_larva>0 
        larva!(du, u, p, t)
    end

    if u.ind.is_metamorph>0
        metamorph!(du, u, p, t)
    end
end