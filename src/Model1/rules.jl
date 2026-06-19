function global_rules!(m)
    return nothing
end

@inline function determine_S_max_hist(S, S_max_hist)
    return max(S, S_max_hist)
end

@inline function death_by_loss_of_structure(S, S_max_hist, S_rel_crit, h_S, dt)::Bool
    return ((S/S_max_hist) < S_rel_crit) && (rand() > exp(-h_S * dt))
end

@inline function death_by_aging(age, a_max)::Bool
    return age >= a_max
end

@inline function stochastic_death(h, dt)::Bool
    return rand() > exp(-h * dt)
end


@inline function check_reproduction_period(time_since_last_repro, tau_R)::Bool
    return time_since_last_repro >= tau_R 
end

@inline function calc_num_offspring(R, X_emb_int)::Int64
    return trunc(R / X_emb_int)
end

function individual_rules!(a, m; init_u_ind, gen_p_ind)::Nothing

    @unpack glb,ind = a.u
    du = a.du
    p = a.p

    ind.age += m.aux.dt

    # ======================================================== #
    # Determining the current life stage
    # This is handled by callbacks in the pure-ODE version
    # ======================================================== #

    if a.u.ind.X_emb > 0
        a.u.ind.embryo = 1.
        a.u.ind.larva = 0.
        a.u.ind.metamorph = 0.
        a.u.ind.juvenile = 0.
        a.u.ind.adult = 0.
    else
        if a.u.ind.H < p.ind.H_j1
            a.u.ind.embryo = 0.
            a.u.ind.larva = 1.
            a.u.ind.metamorph = 0.
            a.u.ind.juvenile = 0.
            a.u.ind.adult = 0.
        elseif (a.u.ind.H >= p.ind.H_j1) && (a.u.ind.E_mt >= 0)
            a.u.ind.embryo = 0.
            a.u.ind.larva = 0.
            a.u.ind.metamorph = 1.
            a.u.ind.juvenile = 0.
            a.u.ind.adult = 0.
        elseif (a.u.ind.H >= p.ind.H_j1) && (a.u.ind.E_mt <= 0) && (a.u.ind.H < p.ind.H_p)
            a.u.ind.embryo = 0.
            a.u.ind.larva = 0.
            a.u.ind.metamorph = 0.
            a.u.ind.juvenile = 1.
            a.u.ind.adult = 0.
        elseif (a.u.ind.H >= p.ind.H_j1) && (a.u.ind.E_mt <= 0) && (a.u.ind.H >= p.ind.H_p)
            a.u.ind.embryo = 0.
            a.u.ind.larva = 0.
            a.u.ind.metamorph = 0.
            a.u.ind.juvenile = 0.
            a.u.ind.adult = 1.
        else
            error("Did not meet any of the life stage criteria - check criteria and state variables. \n X_emb=$(a.u.ind.X_emb), H=$(a.u.ind.H).")
        end
    end

    # ======================================================== #
    # Mortality
    # ======================================================== #

    # ---- aging
    # TODO: easiest way to incorporate temp-dependency into aging is via "temperature-age"
    #       if we actually fit the model to aging data, we could extend to use proposed aging rules from the DEBkiss  

    if death_by_aging(ind.age, p.ind.aux.a_max)
        ind.aux.cause_of_death = 1.
        glb.aging_mortality += 1
    end


    # ---- starvation (loss of structure)

    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    ind.aux.S_max_hist = determine_S_max_hist(ind.S, ind.aux.S_max_hist)

    if death_by_loss_of_structure(ind.S, ind.aux.S_max_hist, p.ind.aux.S_rel_crit, p.ind.aux.h_S, m.aux.dt)
        ind.aux.cause_of_death = 2.
        glb.starvation_mortality += 1
    end

    is_aquatic = sum([ind.larva, ind.metamorph])>0
    is_terrestric = sum([ind.juvenile, ind.adult])>0

    h_b = begin
        if is_aquatic
            h_b = p.ind.h_b_aq
        elseif is_terrestric
            h_b = p.ind.h_b_ter
        else
            h_b = 0.
        end
    end

    # ---- background moratality
    if stochastic_death(h_b, m.aux.dt)
        ind.aux.cause_of_death = 3.
        glb.background_mortality += 1
    end

    # ======================================================== #
    # Reproduction
    # Assuming a very simple rule: 
    # Constant 
    # ======================================================== #

    # reproduction only occurs if the reproduction period has been exceeded
    if check_reproduction_period(ind.aux.time_since_last_repro, p.ind.aux.tau_R) 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg

        N = calc_num_offspring(ind.R, p.ind.X_emb_int)
        if isnan(N)
            println(ind.R, p.ind.X_emb_int)
        end

        for _ in 1:N
            m.aux.idcount += 1 # increment individual counter
            push!(
                m.individuals, 
                IBM.Individual(
                    m.p, 
                    a.u.glb, 
                    init_u_ind, 
                    gen_p_ind; 
                    cohort = a.u.ind.aux.cohort += 1,
                    id = m.aux.idcount
                )
            )
           
            ind.R -= p.ind.X_emb_int # decrease reproduction buffer
            ind.aux.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
        end
        ind.aux.time_since_last_repro = 0. # reset reproduction period
    # if reproduction period has not been exceeded,
    else
        ind.aux.time_since_last_repro += m.aux.dt # track reproduction period
    end
    
    ind.aux.fX = let fX = 1
        if is_aquatic
            fX = f_X(glb.X_aq, p.glb.V_patch_aq, p.ind.K_X_lrv)
        end
        if is_terrestric
            fX = f_X(glb.X_ter, p.glb.A_patch_ter, p.ind.K_X_juv)
        end
        fX
    end
    
    return nothing
end

  