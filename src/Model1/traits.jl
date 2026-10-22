"""
Calculate maximum structural weight for larvae. 
"""
function get_Smax_lrv(spc)
    return ((spc.kappa_emb * spc.eta_IA * spc.dI_max_lrv)/spc.k_M_emb)^3
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


function get_Ww(Wd::Real, watercontent::Real)
    return Wd/(1 - watercontent)
end


function get_Ww(u::ComponentVector{<:Real}, watercontent::Real)
    
    Wd = u.ind.S + u.ind.E_mt
    Ww = get_ww(Wd, watercontent)

    return Ww
end