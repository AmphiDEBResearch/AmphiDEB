
"""
    calc_S_max_juv(spc::ComponentVector)

Calculates maximum structural mass of juveniles and adults from parameters. 
"""
calc_S_max_juv(spc::ComponentVector) = ((spc.kappa_juv * spc.dI_max_juv * spc.eta_IA) / spc.k_M_juv)^3

"""
    calc_S_max_ad(spc::ComponentVector)

Alias for `calc_S_max_juv`.
"""
calc_S_max_ad(spc::ComponentVector) = calc_S_max_juv(spc)


"""
    calc_a_B(sim::AbstractDataFrame)

Extracts age at birth from simulation output. 
"""
calc_a_B(sim::AbstractDataFrame) = sim[sim.X_emb.<=0,:].t |> minimum


"""
    calc_S_B(sim::AbstractDataFrame)

Extracts structural mass at birth from simulation output. 
"""
calc_S_B(sim::AbstractDataFrame) = sim[sim.t==calc_a_B(sim),:].S[1]


"""
    calc_a_j1(sim::AbstractDataFrame)

Extract age at beginning of metamorphosis, `j1`, from simulation output. <br>
Beginning of metamorphosis is defined as the point when feeding rates start to decline. <br>
This can align with Gosner stage 42 for amphibians, but can also be as early as Gosner 38.
"""
calc_a_j1(sim::AbstractDataFrame) = begin
    dI = vcat(0, diff(sim.I))
    return sim[dI .<0,:].t |> minimum
end

"""
    calc_a_j2(sim::AbstractDataFrame)

Extract age at the end of metamorphosis, based on depletion of metamorphic reserve. <br>
Without providing an additional maturity threshold (and therefore an additional argument to this function), 
this is the only plausible option we currently have.<br>
In that case, the end of metamorphosis is reached when the reserve is depleted.
"""
calc_a_j2(sim::AbstractDataFrame) = begin 
    return sim[(sim.E_mt .<= 0) .&& (sim.X_emb .<= 0),:].t |> minimum
end

"""
    calc_a_j2(sim::AbstractDataFrame, H_j2)


Extract age at the end of metamorphosis, based on a maturity threshold `H_j2`. <br>
"""
calc_a_j2(sim::AbstractDataFrame, H_j2::R) where R >: Real = begin 
    return sim[sim.H .> H_j2,:].t |> minimum
end

"""
    calc_metamorphosis_duration(sim::AbstractDataFrame)

Calculate duration of metamorphosis. 
"""
calc_metamorphosis_duration(sim::AbstractDataFrame) = sim[sim.metamorph .== 1,:].t |> x-> maximum(x) - minimum(x)


"""
    calc_vB_growthrate_lrv(spc::ComponentVector)::Float64

Calculate larval von Bertalanffy growth rate.
"""
function calc_vB_growthrate_lrv(spc::ComponentVector)::Float64
    return (spc.eta_AS_lrv/(3*d_V)) * spc.k_M_lrv
end


"""
    calc_H_eq_lrv(spc::ComponentVector)

Calculate equilibrium maturity for larval stage. 
"""
function calc_H_eq_lrv(spc::AmphiDEB.ComponentVector)
    S_max_lrv =  ((spc.kappa_emb * spc.dI_max_lrv * spc.eta_IA) / spc.k_M_emb)^3
    return ((1 - spc.kappa_emb) * spc.dI_max_lrv * spc.eta_IA * S_max_lrv^(2/3)) / spc.k_J_emb
end


"""
    calc_H_eq_juv(spc::ComponentVector)

Calculate equilibrium maturity for juvenile/adult stages. 
"""
function calc_H_eq_juv(spc::AmphiDEB.ComponentVector)
    S_max_juv =  ((spc.kappa_juv * spc.dI_max_juv * spc.eta_IA) / spc.k_M_juv)^3
    return ((1 - spc.kappa_juv) * spc.dI_max_juv * spc.eta_IA * S_max_juv^(2/3)) / spc.k_J_juv
end
