# only difference between Model2 and Model2b is in metamorphosis  (re-setting of maturity)

# we thus import all derivatives from Model2

import ..AmphiDEB.Model2: food_dynamics_firstorder!
import ..AmphiDEB.Model2: y_T
import ..AmphiDEB.Model2: y_T_kap
import ..AmphiDEB.Model2: f_X
import ..AmphiDEB.Model2: calc_S_max
import ..AmphiDEB.Model2: embryo!
import ..AmphiDEB.Model2: larva!
import ..AmphiDEB.Model2: metamorph!
import ..AmphiDEB.Model2: juvenile!
import ..AmphiDEB.Model2: adult!
import ..AmphiDEB.Model2: sys_embryo!
import ..AmphiDEB.Model2: sys_larva!
import ..AmphiDEB.Model2: sys_metamorph!
import ..AmphiDEB.Model2: sys_juvenile!
import ..AmphiDEB.Model2: sys_adult!
import ..AmphiDEB.Model2: sim_embryo
import ..AmphiDEB.Model2: sim_juvenile
import ..AmphiDEB.Model2: sim_adult
import ..AmphiDEB.Model2: isoutofdomain


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