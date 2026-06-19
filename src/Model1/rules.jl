
@inline function determine_S_max_hist(S, S_max_hist)
    return max(S, S_max_hist)
end

@inline function death_by_loss_of_structure(S, S_max_hist, S_rel_crit, h_S, dt)::Bool
    return ((S/S_max_hist) < S_rel_crit) && (rand() > exp(-h_S * dt))
end

@inline function death_by_aging(age, a_max)::Bool
    return age >= a_max
end

@inline function death_by_GUTS(h, dt)::Bool
    return rand() > exp(-h * dt)
end

@inline function check_reproduction_period(time_since_last_repro, tau_R)::Bool
    return time_since_last_repro >= tau_R 
end

@inline function calc_num_offspring(R, X_emb_int)::Int64
    return trunc(R / X_emb_int)
end

function default_individual_rules!(a, m)::Nothing

    @unpack glb,ind = a.u
    du = a.du
    p = a.p

    ind.age += m.dt

    if a.u.ind.X_emb > 0
        a.u.ind.is_embryo = 1.
        a.u.ind.is_larva = 0.
        a.u.ind.is_metamorph = 0.
        a.u.ind.is_juvenile = 0.
        a.u.ind.is_adult = 0.
    else
        if a.u.ind.H < p.u.ind.H_j1
            a.u.ind.is_embryo = 0.
            a.u.ind.is_larva = 1.
            a.u.ind.is_metamorph = 0.
            a.u.ind.is_juvenile = 0.
            a.u.ind.is_adult = 0.
        elseif (a.u.ind.H >= p.u.ind.H_j1) && (a.u.ind.E_mt >= 0)
            a.u.ind.is_embryo = 0.
            a.u.ind.is_larva = 0.
            a.u.ind.is_metamorph = 1.
            a.u.ind.is_juvenile = 0.
            a.u.ind.is_adult = 0.
        elseif (a.u.ind.H >= p.u.ind.H_j1) && (a.u.ind.E_mt <= 0) && (a.u.ind.H < p.u.ind.H_p)
            a.u.ind.is_embryo = 0.
            a.u.ind.is_larva = 0.
            a.u.ind.is_metamorph = 0.
            a.u.ind.is_juvenile = 1.
            a.u.ind.is_adult = 0.
        elseif (a.u.ind.H >= p.u.ind.H_j1) && (a.u.ind.E_mt <= 0) && (a.u.ind.H < p.u.ind.H_p)


        end
    end

    # death due to aging
    if death_by_aging(ind.age, p.ind.a_max)
        ind.cause_of_death = 1.
        glb.aging_mortality += 1
    end

    # life-stage transitions are part of ODE in amphibian model and omitted here
    
    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    ind.S_max_hist = determine_S_max_hist(ind.S, ind.S_max_hist)

    if death_by_loss_of_structure(ind.S, ind.S_max_hist, p.ind.S_rel_crit, p.ind.h_S, m.dt)
        ind.cause_of_death = 2.
        glb.starvation_mortality += 1
    end

    # mortality caused by GUTS submodule, including background mortality
    if death_by_GUTS(ind.h_z, m.dt)
        ind.cause_of_death = 3.
        glb.GUTS_mortality += 1
    end

    # --- reproduction based on a constant reproduction period

    # reproduction only occurs if the reproduction period has been exceeded
    if check_reproduction_period(u.ind.aux.time_since_last_repro, p.ind.aux.tau_R) 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg

        N = calc_num_offspring(u.ind.R, p.ind.X_emb_int)
        if isnan(N)
            println(u.ind.R, p.ind.X_emb_int)
        end

        for _ in 1:N
            m.aux.idcount += 1 # increment individual counter
            push!(
                m.individuals, 
                Individual(
                    m.p, 
                    a.u.glb, 
                    init_u_ind, 
                    gen_p_ind; 
                    cohort = a.u.ind.aux.cohort += 1,
                    id = m.aux.idcount
                )
            )
           
            u.ind.R -= p.ind.X_emb_int # decrease reproduction buffer
            u.ind.aux.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
        end
        u.ind.aux.time_since_last_repro = 0. # reset reproduction period
    # if reproduction period has not been exceeded,
    else
        u.ind.aux.time_since_last_repro += m.aux.dt # track reproduction period
    end
    
    a.u.ind.aux.fX = f_X(u.glb.X, p.glb.V_patch, p.ind.K_X)

    
    return nothing
end

  