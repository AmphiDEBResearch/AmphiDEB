# traits.jl
# calculation of life-history traits and other implied properties from simulation output and parameters 

"""
Calculate maximum structural weight for larvae. 
"""
function get_Smax_lrv(spc)
    return ((spc.kappa_emb * (1 - spc.gamma) * spc.eta_IA * spc.dI_max_lrv)/spc.k_M_emb)^3
end

"""
Calculate absolute maximum assimilation rate (mg/d) for larvae.
"""
function get_dAmax_lrv(spc)
    S_max = get_Smax_lrv(spc)
    return spc.dI_max_lrv * spc.eta_IA * S_max^(2/3)
end

"""
Calculate maximum metamorphic reserve buffer `E_mt_max`.
"""
function get_Hmax(spc)
    dAmax = get_dAmax_lrv(spc)
    return ((1 - spc.kappa_emb) * dAmax)/spc.k_J_emb
end

"""
Get state at birth as `DataFrameRow`.
"""
function get_birth(sim::AbstractDataFrame) 
    larva = sim[sim.larva .== 1,:] 
    if nrow(larva)>0
        return larva[1,:]
    else
        return nothing
    end
end

"""
Get state at climax as `DataFrameRow`.
"""
function get_climax(sim::AbstractDataFrame) 
    metam = sim[sim.metamorph .== 1,:] 
    if nrow(metam)>0
        return metam[1,:]
    else
        return nothing
    end
end

"""
Get state at froglet emergence as `DataFrameRow`.
"""
function get_emergence(sim::AbstractDataFrame)
    juv = sim[sim.juvenile .== 1,:] 
    if nrow(juv)>0
        return juv[1,:]
    else
        return nothing
    end
end

"""
Get wet weight at birth.
"""
function get_Wwb(sim, p)

    birth = get_birth(sim)

    if isnothing(birth)
        return NaN
    end

    Wd_b = birth.S + birth.E_mt # dry weight at birth
    Ww_b = Wd_b / (1 - p.spc.watercontent_larvae) # wet weight at birth

    return Ww_b
end

function get_Wwj1(sim, p)

    climax = get_climax(sim)

    if isnothing(climax)
        return NaN
    end

    Wd_j = climax.S + climax.E_mt
    Ww_j = Wd_j / (1 - p.spc.watercontent_larvae) 

    return Ww_j
end

function get_Wwj2(sim, p)

    emergence = get_emergence(sim)

    if isnothing(emergence)
        return NaN
    end

    Wd_j = emergence.S + emergence.E_mt 
    Ww_j = Wd_j / (1 - p.spc.watercontent_juveniles) 

    return Ww_j
end

function get_traits_birth(sim, p)

    birth = get_birth(sim)

    if isnothing(birth)
        return (
            Ww_b = NaN, t_b = NaN
        )
    end

    Wd_b = birth.S + birth.E_mt # dry weight at birth
    Ww_b = Wd_b / (1 - p.spc.watercontent_larvae) # wet weight at birth
    t_b = birth.t

    return (Ww_b = Ww_b, t_b = t_b)
end

function get_traits_climax(sim, p)

    climax = get_climax(sim)

    if isnothing(climax)
        return (
            Ww_j1 = NaN, t_j1 = NaN
        )
    end

    Wd_j = climax.S + climax.E_mt
    Ww_j = Wd_j / (1 - p.spc.watercontent_larvae) 
    t_j = climax.t

    return (Ww_j1 = Ww_j, t_j1 = t_j)
end

function get_traits_emergence(sim, p)

    emergence = get_emergence(sim)

    if isnothing(emergence)
        return (
            Ww_j2 = NaN, t_j2 = NaN
        )
    end

    Wd_j = emergence.S + emergence.E_mt
    Ww_j = Wd_j / (1 - p.spc.watercontent_juveniles) 
    t_j = emergence.t

    return (Ww_j2 = Ww_j, t_j2 = t_j)
end
