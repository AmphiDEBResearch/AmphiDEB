using Pkg; Pkg.activate("test/")

using EcotoxSystems, AmphiDEB
using Plots, StatsPlots
using DataFrames

begin
    # FIXME: what happens when gamma is high?
    #    - 
    p = deepcopy(AmphiDEB.defaultparams)
    p.glb.t_max = 365.  
    p.glb.dX_in = 1e10

    p.spc.H_j1 = Inf
    
    sims = DataFrame()

    for gamma in 0.1:0.1:0.8
        p.spc.gamma = gamma
        sim = AmphiDEB.ODE_simulator(p)
        sim[!,:drymass] = sim.S += sim.E_mt
        sim[!,:drymass_scaled] = sim.drymass ./ maximum(sim.drymass)
        sim[!,:t_scaled] = sim.t ./ maximum(sim.t)
        sim[!,:gamma] .= gamma
        append!(sims, sim)
    end
    
    @df sims plot(
        plot(
            :t_scaled, :drymass_scaled, group = :gamma, 
            palette = palette(:viridis, length(unique(sims.gamma))), 
            lw = 1.5, 
            xlabel = "scaled time", 
            ylabel = "scaled mass"
            )
        #plot(:t, vcat(0, diff(:I)))
    )
            
end

