birth_condition(u, t, integrator) = u.ind.X_emb
function birth_affect!(integrator)
    integrator.u.ind.is_embryo = 0.
    integrator.u.ind.is_juvenile = 1.
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
    u.ind.H - integrator.p.ind.H_j # NOTE: for chemical effects, make sure that the effect is applied on H_j here
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
    (u.ind.E_mt <= 0) && (u.ind.H >= p.ind.H_j) && (u.ind.H < p.ind.H_p)
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
    y_M = 1., y_MP = 1.,
    y_A = 1., y_AP = 1.
    )::Nothing

    @unpack dI_max_lrv, eta_IA, eta_AS_emb, eta_SA, kappa_emb, k_M_emb, k_J_emb, b_T, T_ref, T_A. K_X_aq = p.ind
    @unpack T, V_patch_aq = p.glb
    @unpack X_aq = u.glb

    yT = y_T(T_A, T_ref, T)
    fX = f_X(X_aq, V_patch_aq, K_X_aq)
    kappa_T = y_T_kap(kappa_emb, b_T, T_ref, T)

    dI = fX * dI_max_lrv * S^(2/3) * yT
    dA = dI * eta_IA * y_A * y_AP
    dM = S * k_M * y_M * y_MP * y_T
    dJ = H * k_J * y_M * y_MP * y_T

    dS = Base.ifelse(
        kappa_T * dA >= dM, 
        y_G * y_GP * eta_AS_emb * (kappa_T * dA - dM),
        -(dM / eta_SA - kappa_T * dA)
    )

    dE_mt = Base.ifelse(
        (kappa * dA) > dM, 
        eta_AS * y_G * y_G_P * gamma * (kappa * dA - dM),
        -(dM / eta_SA - (1 - gamma) * kappa * dA)/delta_E
    )

    dE_mt_max = dE_mt

    dH = max(0, (1 - kappa) * dA - dJ)

    du.ind.X_emb = 0.
    du.ind.I = dI
    du.ind.A = dA
    du.ind.M = dM
    du.ind.J = dJ
    du.ind.S = dS
    du.ind.E_mt = dE_mt
    du.ind.dE_mt_max = dE_mt_max
    du.ind.R = 0.
    du.ind.H = dH

    return nothing
end


function metamorph!(du, u, p, t)::Nothing

    return nothing
end


function juvenile!(du, u, p, t)::Nothing

    return nothing
end

function adult!(du, u, p, t)::Nothing

    return nothing
end