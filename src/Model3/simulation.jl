# ======================================================== #
# ODE
# ======================================================== #

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

# ======================================================== #
# IBM
# ======================================================== #