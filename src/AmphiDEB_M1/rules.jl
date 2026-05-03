
"""
    AmphiDEB_global_rules!(m)

Global rule-based portion of the AmphiDEB IBM. 
"""
function AmphiDEB_global_rules!(m)
    m.u.glb.N = length(m.individuals) # tracking population size
    m.u.glb.X = max.(0, m.u.glb.X) # HOTFIX : negative resource abundances can cause chaos
    
    # record counts per life stage
    m.u.glb.N_emb = sum(map(x -> isapprox(1, x.u.ind.embryo, atol = 0.01), m.individuals))
    m.u.glb.N_lrv = sum(map(x -> isapprox(1, x.u.ind.larva, atol = 0.01), m.individuals))
    m.u.glb.N_mt = sum(map(x -> isapprox(1, x.u.ind.metamorph, atol = 0.01), m.individuals))
    m.u.glb.N_juv = sum(map(x -> isapprox(1, x.u.ind.juvenile, atol = 0.01), m.individuals))
    m.u.glb.N_ad = sum(map(x -> isapprox(1, x.u.ind.adult, atol = 0.01), m.individuals))
    
    m.u.glb.N_emb = sum(isapprox(1, x.u.ind.embryo; atol=0.01) for x in m.individuals)
    m.u.glb.N_lrv = sum(isapprox(1, x.u.ind.larva;  atol=0.01) for x in m.individuals)
    m.u.glb.N_mt  = sum(isapprox(1, x.u.ind.metamorph; atol=0.01) for x in m.individuals)
    m.u.glb.N_juv = sum(isapprox(1, x.u.ind.juvenile; atol=0.01) for x in m.individuals)
    m.u.glb.N_ad  = sum(isapprox(1, x.u.ind.adult;    atol=0.01) for x in m.individuals)
end

@inline function determine_S_max_hist(
    S::Float64,
    S_max_hist::Float64
    )::Float64

    return max(S, S_max_hist)

end

@inline function death_by_loss_of_structure(
    S::Float64,
    S_max_hist::Float64,
    S_rel_crit::Float64,
    h_S::Float64,
    dt::Float64
    )::Bool

    return ((S/S_max_hist) < S_rel_crit) && (rand() > exp(-h_S * dt))

end

@inline function death_by_aging(
    age::Float64,
    a_max::Float64
    )::Bool

    return age >= a_max

end

@inline function death_by_GUTS(
    h::Float64, 
    dt::Float64
    )::Bool

    return rand() > exp(-h * dt)

end

@inline function check_reproduction_period(time_since_last_repro::Float64, tau_R::Float64)::Bool
    return time_since_last_repro >= tau_R 
end

@inline function calc_num_offspring(R::Float64, X_emb_int::Float64)::Int64
    return trunc(R / X_emb_int)
end

function default_individual_rules!(
<<<<<<< HEAD:src/AmphiDEB_M1/rules.jl
    a::EcotoxSystems.AbstractIndividual, 
=======
    a::EcotoxSystems.Individual, 
>>>>>>> 92b859b (restructured embryo and larva derivatives):src/individual_rules.jl
    m::EcotoxSystems.IndividualBasedModel
    )::Nothing

    @unpack glb,ind = a.u
    du = a.du
    p = a.p

    ind.age += m.dt

    # death due to aging
    if death_by_aging(ind[:age], p.ind[:a_max])
        ind.cause_of_death = 1.
        glb.aging_mortality += 1
    end

    # life-stage transitions are part of ODE in amphibian model and omitted here
    
    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    ind.S_max_hist = determine_S_max_hist(ind[:S], ind[:S_max_hist])

    if death_by_loss_of_structure(ind[:S], ind[:S_max_hist], p[:ind][:S_rel_crit], p.ind[:h_S], m.dt)
        ind.cause_of_death = 2.
        glb.starvation_mortality += 1
    end

    # mortality caused by GUTS submodule, including background mortality
    if death_by_GUTS(ind[:h_z], m.dt)
        ind.cause_of_death = 3.
        glb.GUTS_mortality += 1
    end

    # reproduction, assuming a constant reproduction period

    # reproduction only occurs if the reproduction period has been exceeded
    if check_reproduction_period(ind[:time_since_last_repro], p.ind[:tau_R]) 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg
        for _ in 1:calc_num_offspring(ind[:R], p[:ind][:X_emb_int])
            m.idcount += 1 # increment individual counter
            push!(m.individuals, EcotoxSystems.DEBIndividual( # create new individual and push to individuals vector
                m.p, 
                m.u.glb; 
                id = m.idcount, 
                cohort = Int(ind.cohort + 1),
                individual_ode! = a.individual_ode!,
                individual_rules! = a.individual_rules!,
                init_individual_statevars = a.init_individual_statevars,
                gen_ind_params = a.generate_individual_params,
                )
            )
            ind.R -= p[:ind][:X_emb_int] # decrease reproduction buffer
            ind.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
        end
        ind.time_since_last_repro = 0. # reset reproduction period
    # if reproduction period has not been exceeded,
    else
        ind.time_since_last_repro += m.dt # track reproduction period
    end

    return nothing
end