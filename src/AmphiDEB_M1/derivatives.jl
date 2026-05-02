# derivatives_M1.jl

const PMOAS = [
    "G", # decrease in growth efficiency
    "M", # increase in maintenance costs
    "A", # decrease in assimilation efficiency
    "R", # decrease in reproduction efficiency
    "Hneg", # decrease in threshold for metamorphosis
    "Hpos", # increase in threshold for metamorphosis
    "KAPneg" # decrease in kappa - faster development, less growth
    ]

"""
    M1_complete_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

ODE system with log-logistic toxicodynamics. 
"""
function M1_complete_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

    AmphiDEB_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    M1_individual_ODE_with_loglogistic_TD!(du, u, p, t)

    return nothing
end

# linear toxicodynamics (classic "DEBtox" stress function)

"""
    M1_complete_ODE_with_linear_TD!(du, u, p, t)::Nothing

ODE system with linear toxicodynamics. 
"""
function M1_complete_ODE_with_linear_TD!(du, u, p, t)::Nothing

    AmphiDEB_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    M1_individual_ODE_with_linear_TD!(du, u, p, t)

    return nothing
end

### individual-level ODE configurations

function M1_individual_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

    pathogen_dynamics!(du, u, p, t)::Nothing 
    TKTD_LL2!(du, u, p, t) # TKTD with mixtures assuming IA
    Pathogen_Infection!(du, u, p, t) # infection, release of zoospores and relative response to sporangia density
    M1_DEB!(du, u, p, t) # Amphibian DEB model

    return nothing
end

function M1_individual_ODE_with_linear_TD!(du, u, p, t)::Nothing

    pathogen_dynamics!(du, u, p, t)::Nothing 
    TKTD_linear!(du, u, p, t) # TKTD with mixtures assuming IA
    Pathogen_Infection!(du, u, p, t) # infection, release of zoospores and relative response to sporangia density
    M1_DEB!(du, u, p, t) # Amphibian DEB mode

    return nothing

end

""" 
    Amphibian_DEB!(du, u, p, t)::Nothing

Definition of the amphibian DEB model M1 (without TKTD).

The default amphibian DEB model assumes a metamorphic reserve compartment,
which is accumulated during larval development, and depleted during metamorphosis. 
"""
function M1_DEB!(du, u, p, t)::Nothing

    determine_life_stage!(du, u, p, t)
    Arrhenius!(du, u, p, t)
    eta_AS, kappa = life_stage_and_plasticity_effects!(du, u, p, t)

    M1_ingestion!(du, u, p, t)
    maintenance!(du, u, p, t)
    growth!(du, u, p, t, eta_AS, kappa)
    maturation!(du, u, p, t, kappa)
    M1_metamorphic_reserve!(du, u, p, t, eta_AS, kappa)

    reproduction!(du, u, p, t, kappa)

    return nothing
end

## callbacks (events)

# adding pathogen spores at a fixed simulation time
condition_inoculation(u, t, integrator) = integrator.p.glb.pathogen_inoculation_time - t
function effect_inoculation!(integrator) 
    integrator.u.glb.P_Z = integrator.p.glb.pathogen_inoculation_dose
end

# medium renewal, result in removal of spores
condition_renewal(u, t, integrator) = prod(integrator.p.glb.medium_renewals  .- t)  
function effect_renewal!(integrator)
    integrator.u.glb.P_Z = 0.
end

# function to generate define all callbacks in a single CallbackSet
"""
    AmphODE_callbacks()

This function generates a default set of callbacks, i.e. discrete events. 

To simulate custom events, one can write a modified version of this function and provide its output as keyword argument `callbacks` 
to `ODE_simulator`. 
"""
function AmphODE_callbacks()

    cb_inoculation = ContinuousCallback(condition_inoculation, effect_inoculation!)
    cb_renewal = ContinuousCallback(condition_renewal, effect_renewal!)

    return CallbackSet(
        cb_inoculation, 
        cb_renewal
        )

end

## derivative functions

### global derivatives

function AmphiDEB_global!(du, u, p, t)::Nothing

    @inbounds begin
        for i in eachindex(u.glb[:X])
            du.glb.X[i] = dX(p.glb[:dX_in][i], p.glb[:k_V][i], u.glb[:X][i])
        end        
    end

    #u.ind.T = p.glb.tempfun(t, p.glb.temp...) # TODO: temporally variable temperature could be added here, e.g. applying temperature_sinusoidal

    return nothing
end


"""
    temperature_sinusoidal(t::Real, T_max::Real, T_min::Real, t_peak ::Real)::Real

Calculate seasonal fluctuations in temperature from a sinusoidal function.
Currently not included in the default model.
"""
@inline function temperature_sinusoidal(t::Real, T_mean::Real, T_amp::Real, T_phi::Real)::Real

    return T_mean + T_amp * sin*(2π /365 * t + T_phi)

end

"""
    dX(
        dX_in::Real, 
        k_V::Real, 
        X::Real
        )::Real

Change in external resource abundace `X`. 

- `dX_in`: Resource input rate [m/t]
- `k_V`: Resource dilution or mortality rate [1/d]
- `X`: Resource abundance [m]
"""
@inline function dX(
    dX_in::Real, 
    k_V::Real, 
    X::Real
    )::Real

    return dX_in - k_V * X

end

"""
    function dP_Z(
        mu::Real,
        P_Z::Real
        )::Real

Derivative of the zoospore abundance `P_Z`. 

Note that this funcpathogen_dynamics!(du, u, p, t)::Nothing 

endtion does not take changes in zoopsore abundance due to encystem into account. 
This is handled in `Pathogen_Infection!`.

- `mu`: Zoospore background mortality rate
- `P`_Z`: Current zoospore abundance in the environment
"""
@inline function dP_Z(
    mu::Real,
    P_Z::Real
    )::Real

    return -mu * P_Z

end

function Pathogen_growth!(du, u, p, t)::Nothing
    
    du.glb.P_Z = dP_Z(p.pth[:mu], u.glb[:P_Z]) # change in zoospores in environment
    
    return nothing
end

## individual-level derivatives

"""
    function dP_S(
        v0::Real,
        gamma::Real,
        P_Z::Real,
        eta::Real,
        f::Real,
        P_S::Real,
        sigma0::Real,
        sigma1::Real,
        Chi::Real
        )::Real

Derivative of the individual-specific sporangia abundace `P_S`.
See [Drawert et al (2017)](https://royalsocietypublishing.org/doi/full/10.1098/rsif.2017.0480) for details on parameters and model formulation.
"""
@inline function dP_S(
    v0::Real,
    gamma::Real,
    P_Z::Real,
    eta::Real,
    f::Real,
    P_S::Real,
    sigma0::Real,
    sigma1::Real,
    Chi::Real
    )::Real

    return v0 * gamma * P_Z + v0 * eta * f * P_S - (sigma0 + sigma1 * Chi * P_S) * P_S

end

function pathogen_dynamics!(du, u, p, t)::Nothing 

    # growth and killing of sporangia 
    du.ind.P_S = dP_S(p.pth[:v0], p.pth[:gamma], u.glb[:P_Z], p.pth[:eta], p.pth[:f], u.ind[:P_S], p.pth[:sigma0], p.pth[:sigma1], p.ind[:Chi])

    # feedback with zoospore population 
    du.glb.P_Z += p.pth[:eta] * (1-p.pth[:f]) * u.ind[:P_S] # release of spores
    du.glb.P_Z -= p.pth[:gamma] * u.glb[:P_Z] # encystment

    return nothing

end

function pathogen_response!(du, u, p, t)::Nothing
    
    # relative response to pathogen 
    @inbounds begin
        E_P = p.ind[:E_P]
        B_P = p.ind[:B_P]
        for j in eachindex(u.ind[:y_jP])
            if j != 2
                u.ind.y_jP[j] = LL2(u.ind[:P_S], E_P[j], B_P[j])
            else
                u.ind.y_jP[j] = LL2pos(u.ind[:P_S], E_P[j], B_P[j])
            end
        end    
    end

   
    return nothing
end

function Pathogen_Infection!(du, u, p, t)::Nothing

    pathogen_dynamics!(du, u, p, t)
    pathogen_response!(du, u, p, t)

    return nothing
end

#### TKTD functions

@inline function LL2(x::Real, p::NTuple{2,Float64})::Real
    return (1 / (1 + Complex(x / p[1]) ^ p[2])).re
end

@inline function LL2pos(x::Real, p::NTuple{2,Float64})::Real
    return 1 - log((1 / (1 + Complex(x / p[1]) ^ p[2])).re)
end

@inline function NEC2pos(x::Real, p::NTuple{2,Float64})::Real
    return 1 + (p[2] * max(0, x - p[1]))
end

@inline function NEC2neg(x::Real, p::NTuple{2,Float64})::Real
    return min(1, 1 - (p[2] * max(0, x - p[1])))
end

@inline LL2(x::Real, p1::Real, p2::Real)::Real = LL2(x, (p1, p2))
@inline LL2pos(x::Real, p1::Real, p2::Real)::Real = LL2pos(x, (p1, p2))
@inline LL2GUTS(x::Real, p1::Real, p2::Real)::Real = -log(LL2(x, (p1, p2)))


@inline NEC2neg(x::Real, p1::Real, p2::Real)::Real = NEC2neg(x, (p1, p2))
@inline NEC2pos(x::Real, p1::Real, p2::Real)::Real = NEC2pos(x, (p1, p2))


"""
    TK_aquatic(
        larva::Real,
        k_D::Real, 
        C_W::Real,
        D::Real,
        fb_u::Real, 
        fb_e::Real, 
        fb_g::Real,
        fb_R::Real
        )::Real

TK model for aquatic exposure of amphibians. 
Aquatic exposure is 0 for non-larvae.

## Arguments 

- `larva`: Indicator of whether current life stage is larval
- `k_D`: Dominant rate constant 
- `C_W`: Aquatic exposure concentration
- `D`: Scaled damage
- `chi_G`: TK feedback with growth.

## References 

Jager, T. (2020). Revisiting simplified DEBtox models for analysing ecotoxicity data. Ecological Modelling, 416, 108904.
"""
@inline function TK_aquatic(
    larva::Real,
    k_D::Real, 
    C_W::Real,
    D::Real,
    chi_G::Real
    )::Real

    return k_D * ((larva * C_W) - D) - chi_G * D

end

function calc_chi_u(
    fb_u::Real, 
    S_max::Real, 
    S::Real, 
    )::Real

    return fb_u * S_max^(1/3)/S

end

function calc_chi_e(
    fb_e::Real, 
    S_max::Real, 
    S::Real
    )::Real

    return fb_e * S_max^(1/3)/S

end

function calc_chi_G(
    fb_G::Real, 
    dS::Real,
    dE_mt::Real, 
    S::Real,
    E_mt::Real,
    )::Real


    return fb_G * (dS+dE_mt)/(S+E_mt)

end


"""
    TKTD_LL2!(du, u, p, t)::Nothing

TKTD module with log-logistic dose-response. 
"""
@inline function TKTD_LL2!(du, u, p, t)::Nothing
    
    @unpack glb, ind = u

    ind.y_j .= 1.0 # reset relative responses 
    ind.h_z = p.ind[:h_b] # reset GUTS-SD hazard rate to background mortality

    # calculate feedbacks
    ind.chi_G = calc_chi_G(
        p.ind[:fb_G], 
        du.ind[:S], 
        du.ind[:E_mt], 
        u.ind[:S], 
        u.ind[:E_mt]
        )

    @inbounds begin
        for z in eachindex(glb.C_W) # for every chemical
            for j in eachindex(ind.y_j) # for every PMoA
                # calculate change in damage
                du.ind.D_j[z,j] = TK_aquatic(
                    ind[:larva], 
                    p.ind.KD[z,j], 
                    glb.C_W[z], 
                    ind.D_j[z,j],
                    ind[:chi_G]
                    ) 
                # update relative response with respect to PMoA j
                # PMoAs with decreasing response
                if !(j in [2,6]) 
                    ind.y_j[j] *= LL2(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
                # PMoAs with increasing response
                else
                    ind.y_j[j] *= LL2pos(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
                end
            end
            # calculate change in damage for lethal effects
            du.ind.D_h[z] = TK_aquatic(ind[:larva], p.ind[:KD_h][z], glb[:C_W][z], ind[:D_h][z], ind[:chi_G]) 
            # update hazard rate
            ind.h_z += LL2GUTS(ind.D_h[z], p.ind.E_h[z], p.ind.B_h[z])
        end
    end

    du.ind.S_z = -ind.h_z * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

"""
    TKTD_linear!(du, u, p, t)::Nothing

TKTD module with linear dose-response.
"""
@inline function TKTD_linear!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    ind.y_j .= 1.0 # reset relative responses 
    ind.h_z = p.ind[:h_b] # reset GUTS-SD hazard rate to background mortality

    for z in eachindex(glb.C_W) # for every chemical
        for j in eachindex(ind.y_j) # for every PMoA
            # calculate change in damage
            du.ind.D_j[z,j] = TK_aquatic(ind[:larva], p.ind.KD[z,j], glb.C_W[z], ind.D_j[z,j]) 
            # update relative response with respect to PMoA j
            # PMoAs with decreasing response
            if j != 2 
                ind.y_j[j] *= NEC2neg(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
            # PMoAs with increasing response
            else
                ind.y_j[j] *= NEC2pos(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
            end
        end
        # calculate change in damage for lethal effects
        du.ind.D_h[z] = TK_aquatic(ind[:larva], p.ind[:KD_h][z], glb[:C_W][z], ind[:D_h][z]) 
        # update hazard rate
        ind.h_z += LL2GUTS(ind.D_h[z], p.ind.E_h[z], p.ind.B_h[z])
    end

    du.ind.S_z = -ind.h_z * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

"""
    is_embryo(
        X_emb::Real;
        beta::Real = 1e3
        )::Real

Determine whether an individual is in embryonic state, returning life stage indicator as float. 
Following rules of the standard DEBkiss model, embryonic state is maintained until the embryonic buffer is emptied (`X_emb_int`=0).
"""
@inline function is_embryo(
    X_emb::Real;
    beta::Real = 1e3
    )::Real

    return sig(X_emb, 0., 0., 1., beta = beta)

end

"""
    is_larva(
        X_emb::Real,
        H::Real,
        H_j1::Real,
        y_H_neg::Real, 
        y_H_pos::Real;
        beta1::Real = 1e3,
        beta2::Real = 1e7
        )::Real

Determine wheter an individual is in larval state, returning life stage indicator as float. 
Larvae are all non-embryos whose maturity is below the threshold for metamorphosis `H_j1`. 
The timing of `H_j1` will often be aligned with Gosner stage 42 for anurans, but might ocurr earlier. 
The decisive factor is the decline in feeding rates.
"""
@inline function is_larva(
    X_emb::Real,
    H::Real,
    H_j1::Real,
    y_H_neg::Real, 
    y_H_pos::Real;
    beta1::Real = 1e3,
    beta2::Real = 1e7
    )::Real

    return sig(X_emb, 0., 1., 0., beta = beta1) * sig(H, H_j1 * y_H_neg * y_H_pos, 1., 0., beta = beta2)

end

"""
    is_metamorph(
        H::Real, 
        H_j1::Real, 
        y_H_neg::Real,
        y_H_pos::Real,
        E_mt::Real;
        beta1::Real = 1e7, 
        beta2::Real = 1e3, 
        )::Real

Determine wheter an individual is in metamorph state, returning life stage indicator as float.
Metamorphs have a maturity level above ``H_j1`, and their metamorphic reserve `E_mt` is not emptied yet.
"""
@inline function is_metamorph(
    H::Real, 
    H_j1::Real, 
    y_H_neg::Real,
    y_H_pos::Real,
    E_mt::Real;
    beta1::Real = 1e7, 
    beta2::Real = 1e3, 
    )::Real

    return sig(H, H_j1 * y_H_neg * y_H_pos, 0., 1., beta = beta1) * sig(E_mt, 0., 0., 1., beta = beta2) 

end

"""
    is_juvenile(
        H::Real,
        H_j1::Real,
        y_H_neg::Real,
        y_H_pos::Real,
        E_mt::Real,
        H_p::Real;
        beta1 = 1e3,
        beta2 = 1e3,
        beta3 = 1e3
        )::Real

Determine wheter an individual is in juvenile state, returning life stage indicator as float.
Juveniles are individuals with a maturity level between `H_j1` and `H_p` whose metamorphic reserve `E_mt `has been emptied.
"""
@inline function is_juvenile(
    H::Real,
    H_j1::Real,
    y_H_neg::Real,
    y_H_pos::Real,
    E_mt::Real,
    H_p::Real;
    beta1 = 1e3,
    beta2 = 1e3,
    beta3 = 1e3
    )::Real

    return sig(H, H_j1 * y_H_neg * y_H_pos, 0., 1., beta = beta1) * sig(E_mt, 0., 1., 0., beta = beta2) * sig(H, H_p, 1., 0., beta = beta3)

end

"""
    is_adult(
        H::Real,
        H_p::Real;
        beta = 1e3
        )::Real

Determine wheter an individual is in adult state, returning life stage indicator as float.
Adults are individuals with a maturity level above the puberty threshold `H_p`.
"""
@inline function is_adult(
    H::Real,
    H_p::Real;
    beta = 1e3
    )::Real

    return sig(H, H_p, 0., 1., beta = beta)

end


"""
    determine_life_stage!(du, u, p, t)::Nothing


Determine an individual's current life stages by calculating all life stage indicators. 
To avoid numerical issues during ODE solve, all life stage transitions are approximated by sigmoid functions.
"""
function determine_life_stage!(du, u, p, t)::Nothing

    @inbounds begin

        y_j = u.ind[:y_j]

        u.ind.embryo = is_embryo(u.ind[:X_emb])
        u.ind.larva = is_larva(u.ind[:X_emb], u.ind[:H], p.ind[:H_j1], y_j[5], y_j[6])
        u.ind.metamorph = is_metamorph(u.ind[:H], p.ind[:H_j1], y_j[5], y_j[6], u.ind[:E_mt]) 
        u.ind.juvenile = is_juvenile(u.ind[:H], p.ind[:H_j1], y_j[5], y_j[6], u.ind[:E_mt], p.ind[:H_p]) 
        u.ind.adult = is_adult(u.ind[:H], p.ind[:H_p]) 
        
    end

    return nothing
end

"""
    calc_eta_AS(
        embryo::Real,
        larva::Real,
        metamorph::Real,
        eta_AS_emb::Real,
        juvenile::Real,
        adult::Real,
        eta_AS_juv::Real
        )::Real

Individual's current value of the growth efficiency `eta_AS`, 
based on life stage indicators and life stage-specific values of `eta_AS`.
"""
@inline function calc_eta_AS(
    embryo::Real,
    larva::Real,
    metamorph::Real,
    eta_AS_emb::Real,
    juvenile::Real,
    adult::Real,
    eta_AS_juv::Real
    )::Real

    return (embryo + larva + metamorph) * eta_AS_emb + (juvenile + adult) * eta_AS_juv

end


"""
    calc_kappa(
        embryo::Real,
        larva::Real,
        metamorph::Real,
        kappa_emb::Real,
        juvenile::Real,
        adult::Real,
        kappa_juv::Real,
        b_T::Real, 
        T_ref::Real,
        y_K::Real
        )::Real

Individual's current value of the growth efficiency `kappa`, 
based on life stage indicators and life stage-specific values of `kappa`.

This includes temperature effects on κ as in Romoli et al. (2024).
"""
@inline function calc_kappa(
    embryo::Real,
    larva::Real,
    metamorph::Real,
    kappa_emb::Real,
    juvenile::Real,
    adult::Real,
    kappa_juv::Real,
    b_T::Real, 
    T_ref::Real,
    T::Real,
    y_K::Real
    )::Real

    kappa = (embryo + larva + metamorph) * kappa_emb + (juvenile + adult) * kappa_juv

    return 1/(1 + (((1-kappa)/kappa) * exp(-b_T * ((T_ref - T)/T_ref))))*y_K
end

@inline function calc_S_max(
    dI_max::Real, 
    eta_IA::Real, 
    kappa::Real, 
    k_M::Real
    )::Real

    return ((dI_max*eta_IA*kappa)/k_M)^3

end

function calc_S_max(
    embryo::Real, 
    larva::Real, 
    metamorph::Real, 
    juvenile::Real, 
    adult::Real, 
    dI_max_emb::Real, 
    dI_max_lrv::Real,
    dI_max_juv::Real, 
    eta_IA::Real, 
    kappa::Real, 
    k_M_emb::Real, 
    k_M_juv::Real
    )::Real

    return (embryo) * calc_S_max(dI_max_emb, eta_IA, kappa, k_M_emb) +
           (larva + metamorph) * calc_S_max(dI_max_lrv, eta_IA, kappa, k_M_emb) + 
           (juvenile + adult) * calc_S_max(dI_max_juv, eta_IA, kappa, k_M_juv)

end


"""
    life_stage_and_plasticity_effects!(du, u, p, t)::Tuple{Float64,Float64}

Life-stage specificity of parameters., as well as plastic responses to environmental factors (e.g. effect of temperature on κ).

Updates the value of `u.ind[:S_max]`, returns life stage-specific `kappa` and `eta_AS`.
"""
function life_stage_and_plasticity_effects!(du, u, p, t)::Tuple{Float64,Float64}

    @inbounds begin

        y_j = u.ind.y_j

        eta_AS = calc_eta_AS(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:eta_AS_emb], u.ind[:juvenile], u.ind[:adult], p.ind[:eta_AS_juv])
        kappa = calc_kappa(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:kappa_emb], u.ind[:juvenile], u.ind[:adult], p.ind[:kappa_juv], p.ind[:b_T], p.ind[:T_ref], p.glb[:T], y_j[7])
        u.ind.S_max = calc_S_max(
            u.ind[:embryo],
            u.ind[:larva],
            u.ind[:metamorph],
            u.ind[:juvenile],
            u.ind[:adult],
            p.ind[:dI_max_emb],
            p.ind[:dI_max_lrv],
            p.ind[:dI_max_juv],
            p.ind[:eta_IA],
            kappa, 
            p.ind[:k_M_emb],
            p.ind[:k_M_juv]
        )

        return eta_AS, kappa
    end
end



"""
    y_T(
        T_A::Real,
        T_ref::Real,
        T::Real
        )::Real

Temperature correction coefficient `y_T` according to Arrhenius equation.
"""
@inline function y_T(
    T_A::Real,
    T_ref::Real,
    T::Real
    )::Real

    return exp((T_A / T_ref) - (T_A / T))

end

function Arrhenius!(du, u, p, t)::Nothing
    
    u.ind.y_T = y_T(p.ind[:T_A], p.ind[:T_ref], p.glb[:T])

    return nothing

end


"""
    f_X(
        X::Real,
        V_patch::Real,
        K_X::Real
        )::Real

Scaled functional response `f_X` based on Holling Type II functional response.
"""
@inline function f_X(
    X::Real,
    V_patch::Real,
    K_X::Real
    )::Real

    return (X / V_patch) / ((X / V_patch) + K_X)

end

"""
    f_X(
        larva::Real,
        metamorph::Real, 
        juvenile::Real, 
        adult::Real,
        X::Vector{Float64},
        V_patch::Union{Vector{Real},Vector{Float64}}, 
        K_X_lrv::Real,
        K_X_juv::Real
        )::Real

Scaled functional response with account for life stage-specific half saturation cosntant `K_X`.
"""
@inline function f_X(
    larva::Real,
    metamorph::Real, 
    juvenile::Real, 
    adult::Real,
    X::Vector{Float64},
    V_patch::Union{Vector{Real},Vector{Float64}}, 
    K_X_lrv::Real,
    K_X_juv::Real
    )::Real


    return ((larva + metamorph) * f_X(X[1], V_patch[1], K_X_lrv)) + ((juvenile + adult) * f_X(X[2], V_patch[2], K_X_juv))

end

"""
    calc_dI_emb(
        embryo::Real,
        S::Real,
        dI_max_emb::Real,
        y_T::Real
        )::Real

Embryonic buffer uptake rate for embryos.
"""
@inline function calc_dI_emb(
    embryo::Real,
    S::Real,
    dI_max_emb::Real,
    y_T::Real
    )::Real

    return embryo * (Complex(S)^(2/3)).re * dI_max_emb * y_T
    
end

"""
    calc_dI_mt(
        metamorph::Real,
        f_X::Real,
        dI_max_lrv::Real,
        E_mt::Real,
        E_mt_max::Real,
        S::Real,
        y_T::Real
        )::Real

Residual ingestion rate for metamorphs (if any). 
"""
@inline function calc_dI_mt(
    metamorph::Real,
    f_X::Real,
    dI_max_lrv::Real,
    E_mt::Real,
    E_mt_max::Real,
    S::Real,
    y_T::Real
    )::Real

    return metamorph * f_X * dI_max_lrv * (E_mt / E_mt_max) * (Complex(S)^(2/3)).re * y_T

end

"""
    calc_dI_lrv(
        larva::Real,
        f_X::Real,
        dI_max_lrv::Real,
        S::Real,
        y_T::Real,
        )::Real


Food ingestion rate for larvae. 
"""
@inline function calc_dI_lrv(
    larva::Real,
    f_X::Real,
    dI_max_lrv::Real,
    S::Real,
    y_T::Real,
    )::Real

    return larva * f_X * dI_max_lrv * (Complex(S)^(2/3)).re * y_T

end

"""
    calc_dI_juv(
        juvenile::Real,
        adult::Real,
        f_X::Real,
        dI_max_juv::Real,
        S::Real,
        y_T::Real
        )::Real

Food ingestion rate for juveniles and adults.
"""
@inline function calc_dI_juv(
    juvenile::Real,
    adult::Real,
    f_X::Real,
    dI_max_juv::Real,
    S::Real,
    y_T::Real
    )::Real

    return (juvenile + adult) * f_X * dI_max_juv * (Complex(S)^(2/3)).re * y_T

end

"""
    dI(
        dI_emb::Real,
        dI_mt::Real, 
        dI_lrv::Real, 
        dI_juv_ad::Real
        )::Real

Food ingestion rate with account for life stage-specificity.
"""
@inline function dI(
    dI_emb::Real,
    dI_mt::Real, 
    dI_lrv::Real, 
    dI_juv_ad::Real
    )::Real

    return dI_emb + dI_mt + dI_lrv + dI_juv_ad

end

"""
    dA(
        dI::Real, 
        eta_IA::Real, 
        y_A::Real, 
        y_AP::Real
        )::Real

Assimilation flux, following standard DEBkiss model.
"""
@inline function dA(
    dI::Real, 
    eta_IA::Real, 
    y_A::Real, 
    y_AP::Real
    )::Real

    return dI * eta_IA * y_A * y_AP

end

function M1_ingestion!(du, u, p, t)::Nothing

    @inbounds begin

    u.ind.f_X = f_X(u.ind[:larva], u.ind[:metamorph], u.ind[:juvenile], u.ind[:adult], u.glb[:X], p.glb[:V_patch], p.ind[:K_X_lrv], p.ind[:K_X_juv])

    dI_emb = calc_dI_emb(u.ind[:embryo], u.ind[:S], p.ind[:dI_max_emb], u.ind[:y_T]) # ingestion by embryos
    dI_mt = calc_dI_mt(u.ind[:metamorph], u.ind[:f_X], p.ind[:dI_max_lrv], u.ind[:E_mt], u.ind[:E_mt_max], u.ind[:S], u.ind[:y_T]) # residual ingestion flux by metamorphs (which may include tadpoles close to climax)
    dI_lrv = calc_dI_lrv(u.ind[:larva], u.ind[:f_X], p.ind[:dI_max_lrv], u.ind[:S], u.ind[:y_T]) # ingestion by larvae
    dI_juv_ad = calc_dI_juv(u.ind[:juvenile], u.ind[:adult], u.ind[:f_X], p.ind[:dI_max_juv], u.ind[:S], u.ind[:y_T])
    
    du.ind.I = dI(dI_emb, dI_mt, dI_lrv, dI_juv_ad)
    du.ind.X_emb = -dI_emb
    du.glb.X[1] -= (u.ind[:larva] + u.ind[:metamorph]) * du.ind.I
    du.glb.X[2] -= (u.ind[:juvenile] + u.ind[:adult]) * du.ind.I

    # assimilation flux

        y_j = u.ind[:y_j]
        y_jP = u.ind[:y_jP]
        du.ind.A = dA(du.ind[:I], p.ind[:eta_IA], y_j[3], y_jP[3])
    end
    
    return nothing 
end

"""
    determine_k_M(
        embryo::Real,
        larva::Real,
        metamorph::Real,
        k_M_emb::Real,
        delta_k_M_mt::Real,
        juvenile::Real,
        adult::Real,
        k_M_juv::Real
        )::Real

Current value of `k_M`, depending on life stage.
"""
@inline function determine_k_M(
    embryo::Real,
    larva::Real,
    metamorph::Real,
    k_M_emb::Real,
    delta_k_M_mt::Real,
    juvenile::Real,
    adult::Real,
    k_M_juv::Real
    )::Real

    return (embryo + larva) * k_M_emb + metamorph * k_M_emb * delta_k_M_mt + (juvenile + adult) * k_M_juv

end

"""
    determine_k_J(
        embryo::Real,
        larva::Real,
        metamorph::Real,
        k_J_emb::Real,
        juvenile::Real,
        adult::Real,
        k_J_juv::Real
        )::Real

Current value of `k_J`, depending on life stage.
"""
@inline function determine_k_J(
    embryo::Real,
    larva::Real,
    metamorph::Real,
    k_J_emb::Real,
    juvenile::Real,
    adult::Real,
    k_J_juv::Real
    )::Real

    return (embryo + larva + metamorph) * k_J_emb + (juvenile + adult) * k_J_juv

end

"""
    dM(
        S::Real,
        k_M::Real,
        y_M::Real,
        y_MP::Real,
        y_T::Real
        )::Real

Somatic maintenance rate, 
including response to chemicals (`y_M`), pathogens (`y_MP`) and temperature (`y_T`).
"""
@inline function dM(
    S::Real,
    k_M::Real,
    y_M::Real,
    y_MP::Real,
    y_T::Real
    )::Real

    return S * k_M * y_M * y_MP * y_T

end

"""
    dM(
        S::Real,
        k_M::Real,
        y_M::Real,
        y_MP::Real,
        y_T::Real
        )::Real

Maturity maintenance rates, 
including response to chemicals (`y_M`), pathogens (`y_MP`) and temperature (`y_T`).
"""
@inline function dJ(
    H::Real,
    k_J::Real,
    y_M::Real,
    y_MP::Real,
    y_T::Real
    )::Real

    return H * k_J * y_M * y_MP * y_T

end


function maintenance!(du, u, p, t)::Nothing 

    k_M = determine_k_M(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:k_M_emb], p.ind[:delta_k_M_mt], u.ind[:juvenile], u.ind[:adult], p.ind[:k_M_juv])   
    k_J = determine_k_J(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:k_J_emb], u.ind[:juvenile], u.ind[:adult], p.ind[:k_J_juv])   

    @inbounds begin
        y_j =  u.ind[:y_j]
        y_jP = u.ind[:y_jP]

        du.ind.M = dM(u.ind[:S], k_M, y_j[2], y_jP[2], u.ind[:y_T])
        du.ind.J = dJ(u.ind[:H], k_J, y_j[2], y_jP[2], u.ind[:y_T])    
    end
    
    return nothing
end

"""
    dR(
        adult::Real,
        eta_AR::Real,
        y_R::Real,
        y_RP::Real,
        kappa::Real,
        dA::Real,
        dJ::Real
        )::Real

Reproduction rate, including response to chemicals `y_R` and pathogens (`y_RP`).
"""
@inline function dR(
    adult::Real,
    eta_AR::Real,
    y_R::Real,
    y_RP::Real,
    kappa::Real,
    dA::Real,
    dJ::Real
    )::Real

    return adult * clipneg(eta_AR * y_R * y_RP * ((1 - kappa) * dA - dJ))

end

function reproduction!(du, u, p, t, kappa::Real)::Nothing
    
    @inbounds begin
        y_j = u.ind[:y_j]
        y_jP = u.ind[:y_jP]
        du.ind.R = dR(u.ind[:adult], p.ind[:eta_AR] , y_j[4], y_jP[4], kappa, du.ind[:A], du.ind[:J])
    end

    return nothing
end

"""
    calc_dS_emb_juv_ad(
        kappa::Real,
        dA::Real,
        dM::Real,
        eta_SA::Real,
        y_G::Real,
        y_GP::Real,
        eta_AS::Real
        )::Real

Somatic growth rate for embryos, juveniles and adults, 
including response to chemicals (`y_G`) and pathogens (`y_GP`).
"""
@inline function calc_dS_emb_juv_ad(
    kappa::Real,
    dA::Real,
    dM::Real,
    eta_SA::Real,
    y_G::Real,
    y_GP::Real,
    eta_AS::Real
    )::Real

    return sig( 
        kappa * dA, 
        dM, 
        -(dM / eta_SA - kappa * dA), 
        y_G * y_GP * eta_AS * (kappa * dA - dM) 
    )    

end

"""
    calc_dS_lrv(
        kappa::Real,
        dA::Real,
        dM::Real,
        eta_SA::Real,
        gamma::Real,
        eta_AS::Real,
        y_G::Real,
        y_GP::Real
        )::Real

Somatic growth rate for larvae. 
Approach follows DEBkiss model, but with an additional γ-split to divert assimilates towards the metamorphic reserve. 
This is also taken into account in the shrinking equation (`kappa * dA < dM`).
"""
@inline function calc_dS_lrv(
    kappa::Real,
    dA::Real,
    dM::Real,
    eta_SA::Real,
    gamma::Real,
    eta_AS::Real,
    y_G::Real,
    y_GP::Real
    )::Real

    if (kappa * dA) > dM
        return eta_AS * y_G * y_GP * (1 - gamma) * (kappa * dA - dM)
    else
        return -(dM / eta_SA - (1 - gamma) * kappa * dA)
    end
end

"""
    calc_dS_mt(
        metamorph::Real,
        eta_AS::Real,
        y_G::Real,
        y_GP::Real,
        dA::Real
        )::Real

Somatic growth rate for metamorphs. 
"""
@inline function calc_dS_mt(
    metamorph::Real,
    eta_AS::Real,
    y_G::Real,
    y_GP::Real,
    dA::Real
    )::Real

    return metamorph * eta_AS * y_G * y_GP * dA

end


"""
    dS(
        embryo::Real,
        juvenile::Real,
        adult::Real,
        dS_emb_juv_ad::Real,
        larva::Real,
        dS_lrv::Real,
        metamorph::Real,
        dS_mt::Real
        )::Real

Life stage-specific somatic growth rate.
"""
@inline function dS(
    embryo::Real,
    juvenile::Real,
    adult::Real,
    dS_emb_juv_ad::Real,
    larva::Real,
    dS_lrv::Real,
    metamorph::Real,
    dS_mt::Real
    )::Real

    return (embryo + juvenile + adult) * dS_emb_juv_ad + larva * dS_lrv + metamorph * dS_mt

end

function growth!(du, u, p, t, eta_AS::Real, kappa::Real)::Nothing
    
    @inbounds begin

        y_j = u.ind[:y_j]
        y_jP = u.ind[:y_jP]

        dS_emb_juv_ad = calc_dS_emb_juv_ad(kappa, du.ind[:A], du.ind[:M], p.ind[:eta_SA], y_j[1], y_jP[1], eta_AS)
        dS_lrv = calc_dS_lrv(kappa, du.ind[:A], du.ind[:M], p.ind[:eta_SA], p.ind[:gamma], eta_AS, y_j[1], y_jP[1])
        dS_mt = calc_dS_mt(u.ind[:metamorph], eta_AS, y_j[1], y_jP[1], du.ind[:A])
        
        du.ind.S = dS(u.ind[:embryo], u.ind[:juvenile], u.ind[:adult], dS_emb_juv_ad, u.ind[:larva], dS_lrv, u.ind[:metamorph], dS_mt)
    
    end
    
    return nothing
end 

"""
    calc_dE_mt_lrv(
        eta_AS::Real, 
        y_G::Real,
        y_G_P::Real,
        gamma::Real,
        kappa::Real,
        dA::Real,
        dM::Real,
        eta_SA::Real
        delta_E::Real,
        )::Real

Accumulation of metamorphic reserve by larvae.

>[!Note] Why is the growth and shrinking efficiency applied to the reserve? 
>The main reason are the implied effects of chemical stressors when the PMoA is a decrease in growth efficiency. 
>If the growth efficiency only affects structure, the effect on total dry mass is miniscule since reserve accumulation will only be indirectly affected.

## Arguments

- `eta_AS`: growth efficiency
- `y_G`: relative response of growth efficiency to chemical stressors
- `y_G_P`: relative response of growth effficiency to pathogens
- `gamma`: allocation to metamorphic reserve
- `kappa`: allocation to soma incl. metamorphic reserve
- `dA`: assimilation rate
- `dM`: maintenance rate
- `eta_SA`: shrinking efficiency
- `delta_E`: energy density of metamorphic reserve, relative to structure
"""
@inline function calc_dE_mt_lrv(
    eta_AS::Real, 
    y_G::Real,
    y_G_P::Real,
    gamma::Real,
    kappa::Real,
    dA::Real,
    dM::Real,
    eta_SA::Real,
    delta_E::Real
    )::Real

    if (kappa * dA) > dM
        return eta_AS * y_G * y_G_P * gamma * (kappa * dA - dM)
    else
        return -(dM / eta_SA - (1 - gamma) * kappa * dA)/delta_E
    end
end

"""
    calc_dE_mt_mt(
        dH::Real,
        dJ::Real,
        dM::Real
        delta_E::Real
        )::Real

Depletion of metamorphic reserve by metamorphs.
"""
@inline function calc_dE_mt_mt(
    dH::Real,
    dJ::Real,
    dM::Real,
    delta_E::Real
    )::Real

    return -(dH + dJ + dM)/delta_E

end

"""
    dE_mt(
        larva::Real,
        dE_mt_lrv::Real,
        metamorph::Real,
        dE_mt_mt::Real
        )::Real

Life stage-specific dynamics of the metamorphic reserve.
"""
@inline function dE_mt(
    larva::Real,
    dE_mt_lrv::Real,
    metamorph::Real,
    dE_mt_mt::Real
    )::Real

    return larva * dE_mt_lrv + metamorph * dE_mt_mt

end

"""
    dE_mt_max(
        larva::Real,
        dE_mt::Real
        )::Real

Function to record the maximum metamorphic reserve level, reached at the transition from larva to metamorph.
"""
@inline function dE_mt_max(
    larva::Real,
    dE_mt::Real
    )::Real

    return larva * dE_mt

end

function M1_metamorphic_reserve!(du, u, p, t, eta_AS::Real, kappa::Real)::Nothing
    
    @inbounds begin
        y_j = u.ind[:y_j]
        y_jP = u.ind[:y_jP]
        dE_mt_lrv = calc_dE_mt_lrv(p.ind[:eta_AS_emb], y_j[1], y_jP[1], p.ind[:gamma], kappa, du.ind[:A], du.ind[:M], p.ind[:eta_SA], p.ind[:delta_E]) 
        dE_mt_mt = calc_dE_mt_mt(du.ind[:H], du.ind[:J], du.ind[:M], p.ind[:delta_E])
        du.ind.E_mt = dE_mt(u.ind[:larva], dE_mt_lrv, u.ind[:metamorph], dE_mt_mt)
        du.ind.E_mt_max = dE_mt_max(u.ind[:larva], du.ind[:E_mt])
    end

    return nothing
end

"""
    calc_dH(
        kappa::Real,
        dA::Real,
        dJ::Real
        )::Real

Maturation rate following standard DEBkiss model.
"""
@inline function calc_dH(
    kappa::Real,
    dA::Real,
    dJ::Real
    )::Real

    return clipneg((1 - kappa) * dA - dJ)

end

"""
    calc_dH_mt(
        kappa::Real,
        dM::Real,
        dJ::Real
        )::Real

Maturation rate for metamorphs, calculated from the κ-rule and assuming that soma are decoupled from the κ rule during metamorphosis.
"""
@inline function calc_dH_mt(
    kappa::Real,
    dM::Real,
    dJ::Real
    )::Real

    return clipneg((1 - kappa) * dM) / kappa - dJ

end

"""
    dH(
        adult::Real,
        metamorph::Real,
        dH::Real,
        dH_mt::Real
        )::Real

Life stage-specific maturation rate.
"""
@inline function dH(
    adult::Real,
    metamorph::Real,
    dH::Real,
    dH_mt::Real
    )::Real

    return (1 - adult) * (1 - metamorph) * dH + metamorph * dH_mt
end

function maturation!(du, u, p, t, kappa::Real)::Nothing

    # maturation follows kappa-rule for all but metamorphs
    dH_all = calc_dH(kappa, du.ind[:A], du.ind[:J])

    ## for metamorhps, we assume that all but somatic growth is fueled by E_mt, 
    # allowing us to apply the kappa rule to calculate dH from dM and dJ
    dH_mt = calc_dH_mt(kappa, du.ind[:M], du.ind[:J])

    # for adults, there is no maturation => we only need to differentiate between non-adults and metamorph
    du.ind.H = dH(u.ind[:adult], u.ind[:metamorph], dH_all, dH_mt)

    return nothing
end
