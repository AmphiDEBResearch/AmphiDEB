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

using AmphiDEB
using EcotoxSystems


@testset "Pathogen growth without effects" begin 
    global p = deepcopy(AmphiDEB.defaultparams)

    p.glb.t_max = 60
    p.glb.pathogen_inoculation_time = 20.
    p.glb.pathogen_inoculation_dose = 1e3

    p.glb.dX_in = [15., 15.]

    p.spc.Z = Dirac(1.)

    global sim = @replicates AmphiDEB.ODE_simulator(p) 10
        
    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        plot(:t, :S, group = :replicate, ylabel = "S"),
        plot(:t, :E_mt, group = :replicate, ylabel = "E_mt"),
        plot(:t, :P_Z, group = :replicate, ylabel = "P_Z"),
        plot(:t, :P_S, group = :replicate, ylabel = "P_S"),
        xlabel = "t",
        xrotation = 45
        )
    display(plt)

    # check for final median zoospore and sporangia abundance

    sim_end = sim[sim.t .== maximum(sim.t),:]

    @test 500 <= median(sim_end.P_Z) <= 2000
    @test 1 <= median(sim_end.P_S) <= 20
end


@testset "Pathogen growth with effects" begin 

    pmoa_idx = 1
    global p = deepcopy(AmphiDEB.defaultparams)

    p.glb.t_max = 60
    p.glb.pathogen_inoculation_time = 20.
    p.glb.pathogen_inoculation_dose = 1e3

    p.glb.dX_in = [15., 15.]

    p.spc.Z = Dirac(1.)

    p.spc.E_P[pmoa_idx] = 5.
    p.spc.B_P[pmoa_idx] = 2.

    sim = @replicates AmphiDEB.ODE_simulator(p) 10
        
    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        plot(:t, :S, group = :replicate, ylabel = "S"),
        plot(:t, :E_mt, group = :replicate, ylabel = "E_mt"),
        plot(:t, :P_Z, group = :replicate, ylabel = "P_Z"),
        plot(:t, :P_S, group = :replicate, ylabel = "P_S"),
        xlabel = "t", xrotation = 45
        )
    display(plt)

    sim_end = sim[sim.t .== maximum(sim.t),:]

    @test 1 < median(sim_end.S) < 8
end