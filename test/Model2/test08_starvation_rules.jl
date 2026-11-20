using Pkg; Pkg.activate("test")

using EcotoxSystems
using Plots, StatsPlots
using Test

using Revise
using AmphiDEB
 
@testset "Reserve dynamics under starvation" begin
    p.glb.t_max = 365. * 1

    p.glb.dX_in = [.5, 1.] # limit food availability
    p.spc.H_j1 = Inf # ignore metamorphosis
    p.spc.H_p = Inf # ignore adulthood
    
    E_mt_fin = []
    E_mt_rel_fin = []

    for dX_in in [1., .5, .25]
        p.glb.dX_in[1] = dX_in
        sim = AmphiDEB.ODE_simulator(
            p
            )
        
        # verify that E_mt reaches equilibrium
        @test maximum(diff(sim.E_mt)[end-10:end]) <= 1e-4

        push!(E_mt_fin, sim.E_mt[end])
        push!(E_mt_rel_fin, sim.E_mt[end]/sim.S[end])
    end

    # verify that E_mt decreases with decreasing food availability
    @test E_mt_fin == sort(E_mt_fin, rev = true)
    # verify that relative reserve density decreases with decreasing food availability
    @test E_mt_rel_fin == sort(E_mt_rel_fin, rev = true)
            
end




