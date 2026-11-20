using Pkg; Pkg.activate("test")
using Test
using Distributions
using OrdinaryDiffEq

using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

include("testutils.jl")

using DataFrames, DataFramesMeta
using StatsBase

using Revise

@time import AmphiDEB: defaultparams, ODE_simulator, Amphibian_DEB!, AmphiDEB_ODE!
using AmphiDEB
using EcotoxSystems
import EcotoxSystems: DEBODE_global!

import EcotoxSystems: sig
import EcotoxSystems: constrmvec

using Plots 

@testset "Effects of temperature on metamorphosis traits" begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 80.
    p.glb.pathogen_inoculation_time = Inf
    p.glb.dX_in = [20., 20.]
    p.spc.H_p = 50.

    global sims = DataFrame()

    for T in [15., 20., 25.]
        p.glb.T = T + 273.15

        @time sim = AmphiDEB.ODE_simulator(p);

        sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
        sim[!,:dI] = vcat(0, diff(sim.I))
        sim[!,:dI_rel] = sim.dI ./ (sim.S .^(2/3))
        sim[!,:W_tot] = sim.S .+ sim.E_mt 
        sim[!,:T] .= T

        append!(sims, sim)
    end
        
    @df sims plot(:t, :W_tot, group = :T, label = hcat(unique(:T)...), leg = true) |> display

    global metamorphosis_traits = combine(groupby(sims, :T)) do df

        t_G42 = df[isapprox.(1, df.metamorph, atol = 0.99),:].t[1]
        t_G46 = df[isapprox.(1, df.juvenile, atol = 0.99),:].t[1]


        W_G42 = df[isapprox.(1, df.metamorph, atol = 0.99),:].W_tot[1]
        W_G46 = df[isapprox.(1, df.juvenile, atol = 0.99),:].W_tot[1]

        return DataFrame(
            t_G42 = t_G42, 
            t_G46 = t_G46,
            W_G42 = W_G42, 
            W_G46 = W_G46
        )

    end

    @with metamorphosis_traits begin
        # time to G42 should decrease with increasing T
        @test :t_G42 == sort!(:t_G42, rev = true) 
        
        # time to G46 should decrease with increasing T
        @test :t_G46 == sort!(:t_G46, rev = true) 

        # mass at G42 should decrease with increasing T
        @test :W_G42 == sort!(:W_G42, rev = true)
        
        # mass at G46 should decrease with increasing T
        @test :W_G46 == sort!(:W_G46, rev = true) 

    end
end
