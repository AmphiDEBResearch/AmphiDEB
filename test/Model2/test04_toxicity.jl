using Pkg; Pkg.activate("test")
using Test
using Distributions
using OrdinaryDiffEq

using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

include("testutils.jl")

using DataFrames, DataFramesMeta
using StatsBase

using EcotoxSystems
import EcotoxSystems: constrmvec

using Revise

using AmphiDEB
@time import AmphiDEB: defaultparams


@testset "Toxicity with default parameters" begin 
    
    for pmoa_idx in 1:7
    
        global p = deepcopy(defaultparams)
    
        p.glb.t_max = 20.
        #p.glb.pathogen_inoculation_time = Inf

        p.glb.dX_in = [15. 15.]
        p.spc.H_j1 = 0.2
        p.spc.H_p = 40.

        p.spc.KD .= 0.
        p.spc.KD[pmoa_idx] = 1.
        p.spc.E[pmoa_idx] = 2. 
        p.spc.B[pmoa_idx] = 2.

        # check that the syntax for parameter assignment works 
        @test sum(p.spc.KD .> 0) == 1

        @time global sim = exposure(
            AmphiDEB.ODE_simulator, 
            p,
            Matrix(hcat([0.; 1.; 2.;]...)')
            )
            
        sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
        sim[!,:W_tot] = sim.S .+ sim.E_mt 
        
        plt = @df sim plot(
            plot(:t, :S, group = :C_W_1),
            plot(:t, :E_mt, group = :C_W_1),
            plot(:t, :S .+ :E_mt, group = :C_W_1),
            plot(:t, :y_j_1_1, group = :C_W_1, ylim = (-0.01, 1.01)),
            plot(:t, :y_j_1_2, group = :C_W_1, ylim = (0.99, 2)),
            plot(:t, :y_j_1_3, group = :C_W_1, ylim = (-0.01, 1.01)),
            plot(:t, :y_j_1_4, group = :C_W_1, ylim = (-0.01, 1.01)),
            plot(:t, :y_j_1_5, group = :C_W_1, ylim = (-0.01, 1.01)),
            plot(:t, :y_j_1_6, group = :C_W_1, ylim = (0.99, 2)),
            plot(:t, :y_j_1_7, group = :C_W_1, ylim = (-0.01, 1.01)),
            plot(:t, :R, group = :C_W_1),
            plot(:t, :H, group = :C_W_1),
            plot(:t, :larva, group = :C_W_1, linetype = :steppre),
            xrotation = 45, 
            xlabel = "t", ylabel = ["S" "E_mt" "W" "y_G" "y_M" "y_A" "y_R" "y_Hneg" "y_Hpos" "y_κ" "R" "H" "larva"], 
            leg = [:topleft false false false false false false], legendtitle = "C_W", 
            size = (1000,600), bottommargin = 5mm, leftmargin = 5mm
            )

        display(plt)

        # y-column for currently engaged PMoA 
        c = Symbol("y_j_1_$(pmoa_idx)")

        # for every PMoA, verify that the relative responses are calculated correctly
        for idx in 1:7
            # get the associated y-column
            c_idx = Symbol("y_j_1_$(pmoa_idx)")
            # if this is not the currently engaged pmoa
            if !(c_idx==c)
                # verify that y is close to 1 (== no effect)
                @test unique(0.99 .< sim[:,c_idx] .< 1.01) == [true]
            # if this is the currently engaged pmoa
            else
                # the highest concentration is equal to the median effective damage, 
                # so we expect y to bottom out around 0.5 for the given parameters (or max out at 2. for PMoA M)
                if !(pmoa_idx in [2,6])
                    @test isapprox(0.5, minimum(sim[:,c]), atol = 0.1)
                else
                    @test isapprox(1.69, maximum(sim[:,c]), atol = 0.1)
                end
            end
        end

        # verify monotonically increasing effects on R - this should be the case for all PMoAs

        sim_end = @subset(sim, :t .== maximum(:t))
        sort!(sim_end, :R, rev = true)
        @test sort(sim_end.C_W_1) == sim_end.C_W_1
    end
end

import AmphiDEB: PMOAS

# TODO: add proper test for growth-TK-feedback
begin
    p = deepcopy(AmphiDEB.defaultparams)

    p.spc.KD[3] = .1
    p.spc.E[3] = 1. 
    p.spc.B[3] = 2.
    p.spc.fb_G = 1.

    p.spc.H_j1 = Inf

    @time sim_fB = exposure(AmphiDEB.ODE_simulator, p, [0., 0.125, 0.25, 0.5, ]);

    p.spc.fb_G = 0.
    sim_nofB = exposure(AmphiDEB.ODE_simulator, p, [0., 0.125, 0.25, 0.5, ]);
    @df sim_nofB groupedlineplot!(:t, :D_j_1_3, :treatment_id, palette = palette(:default)[1:4], subplot = 1, linestyle = :dash)
    @df sim_nofB groupedlineplot!(:t, :S .+ :E_mt, :treatment_id, palette = palette(:default)[1:4], subplot = 2, linestyle = :dash)


    df = @subset(sim_fB, :treatment_id .== 2)

    df[!,:dW] = vcat(0, diff(df.S))+vcat(0, diff(df.E_mt))
    df[!,:dD] = vcat(0, diff(df.D_j_1_3))

    @df df plot(
        plot(:t, :dW, xlabel = "t", ylabel = "dW"),
        scatter(:dW, :dD, xlabel = "dW", ylabel = "dD")
    )    
end
