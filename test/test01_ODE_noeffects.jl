using Test
using Distributions
using OrdinaryDiffEq

using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

include("testutils.jl")


# ======================================== #
# Submodule setup
# ======================================== #

using AmphiDEB, AmphiDEB.Model1
using EcotoxSystems
using StatsPlots
using DataFrames, DataFramesMeta
using StatsBase

begin # embryo
    p = deepcopy(Model1.params)
    sim_emb, u0lrv, p_ind = Model1.sim_embryo(p)
    @test abs(sim_emb.X_emb[end]) ≈ 0 atol = 1e-3
end

begin # larva
    sim_lrv, u0mt, p_ind = Model1.sim_larva(p_ind, u0lrv)
    @test sim_lrv.H[end] ≈ p.spc.H_j1    
end

<<<<<<< HEAD
#using BenchmarkTools
#AmphiDEB.ODE_simulator(AmphiDEB.defaultparams);
#@benchmark AmphiDEB.ODE_simulator(AmphiDEB.defaultparams)
#
#VSCodeServer.@profview_allocs AmphiDEB.ODE_simulator(AmphiDEB.defaultparams)
=======
begin # metamorph
    sim_mt, u0juv, p_ind = Model1.sim_metamorph(p_ind, u0mt)
    @df sim_mt plot(:t, [:S :E_mt :I], layout = (1,3), marker = true, leg = false)
    @test sim_mt.E_mt[end] ≈ 0. atol = 1e-3
end

begin # juvenile 
    sim_juv, u0ad, p_ind = Model1.sim_juvenile(p_ind, u0juv)
    @df sim_juv plot(:t, :H)
    sim_juv.H[end] ≈ p_ind.ind.H_p
    sim_juv.R[end] == 0.
end
>>>>>>> accb043 (completed new sim_all function and introduced aliases for backwards-compatibility)

using Infiltrator

begin # adult
    p_ind.glb.t_max = 365. * 10 # we use high t_max to check for solver instabilities when kappa*dA ≈ dM
    sim_ad, uend, p_ind = Model1.sim_adult(p_ind, u0ad)
    @df sim_ad plot(:t, [:S :H :R], layout = (1,3), marker = true)
    
    @test unique(diff(sim_ad.H)) == [0.] # maturation should stop
    @test sum(diff(sim_ad.S) .<= 0) == 0 # growth should continue
    @test sim_ad.R[1] == 0. # repro should start at 0
    @test sum(diff(sim_ad.R) .<= 0) == 0 # repro buffer should grow continuously 
end


begin # all
    @time sim = Model1.sim_all(p, saveat = 1)
    @test sum(diff(sim.t) .< 0) == 0 # check that concatenation of time vectors is plausible

    @df sim plot(:t, [:S :H :E_mt :R], layout = (2,2))
end


#using BenchmarkTools
#Model1.sim_all(p); @benchmark Model1.sim_all(p)
#

# ======================================== #
# Previous setup without submodules
# - kept for backwards compatability 
# ======================================== #

using AmphiDEB.OrdinaryDiffEq
sim = AmphiDEB.ODE_simulator(p)

function run_basetest(m; pmod = p->p)

    global p = deepcopy(Model1.params)

    p.glb.t_max = 365*2
    p.glb.pathogen_inoculation_time = Inf
    p.glb.dX_in = [20., 20.]
    p.spc.H_p = 55.

    S_max_anl = AmphiDEB.calc_S_max_juv(p.spc)
    H_max_anl_juv = AmphiDEB.calc_H_eq_juv(p.spc)

    p = pmod(p) # optional modification of parameters

    @time global sim =  (
            p, 
            reltol = 1e-10,
            model = m,
            saveat = 1/24, # we need high-resolution output to verify the solution
            );

    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:dI] = vcat(0, diff(sim.I))
    sim[!,:dI_rel] = sim.dI ./ (sim.S .^(2/3))
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = plot_statevars(
        @subset(sim, :t .< 40), 
        [:S, :H, :E_mt, :R, :X_emb, :J, :dI_rel, :W_tot, :f_X, :embryo, :larva, :metamorph], 
        xrotation = 45
        )
    hline!([p.spc.H_j1], subplot=2, color = :gray, linestyle = :dash)

    display(plt)

    # check final structural mass
    
    @test isapprox(S_max_anl, maximum(sim.S), atol = 5) 

    # check that all life stage indicators max out close to 1
    
    @test 0.99 < maximum(sim.embryo) < 1.01
    @test 0.99 < maximum(sim.larva) < 1.01
    @test 0.99 < maximum(sim.metamorph) < 1.01
    @test 0.99 < maximum(sim.juvenile) < 1.01
    @test 0.99 < maximum(sim.adult) < 1.01

    # check that the sum of life stage indicators is always approximately 1

    sum_indicators = @. sim.embryo + sim.larva + sim.metamorph + sim.juvenile + sim.adult

    @test unique(isapprox.(1, sum_indicators, atol = 1e-3)) == [true]
    
    # verify that the size-specific maintenance rate is approximately constant by back-calculating k_M from the state variables
        
    sim[!,:dM] = vcat(0, diff(sim.M)) ./ vcat(0, diff(sim.t))
    k_M = sim.dM ./ sim.S 

    reldiff_kM = begin
        kmin, kmax = extrema(k_M[50:end]) # we allow for a small "burn in" period: at the beginning, k_M is close to 0 and the relative error will be large 
        kmax / kmin
    end

    @test 0.9 < reldiff_kM < 1.1

    # same for k_J

    sim[!,:dH] = vcat(0, diff(sim.J)) ./ vcat(0, diff(sim.t))
    k_J = sim.dH ./ sim.H 

    reldiff_kJ = begin
        kmin, kmax = extrema(k_J[50:end]) 
        kmax / kmin
    end

    @test 0.9 < reldiff_kJ < 1.1

    # comparing analytically calculated with simulated equilibrium maturity 
    # for this we need to disengage H_p and re-run the simulation

    p.glb.t_max = 365*10
    p.spc.H_p = Inf
    
    @time global sim = AmphiDEB.ODE_simulator(
            p, 
            saveat = 1/24, # we need high-resolution output to verify the solution
            model = m
            );

    @info "
    Analytically caluclated equilibrium maturity: $(round(H_max_anl_juv, sigdigits = 3))
    Simulated maximum maturity: $(round(maximum(sim.H)))
    "

    @test 0.8*H_max_anl_juv <= maximum(sim.H) <= 1.2*H_max_anl_juv 
end

@testset "Default parameters M1" begin
    run_basetest(AmphiDEB.M1_complete_ODE_with_loglogistic_TD!)
end

#
#@testset "Default parameters M2" begin
#    run_basetest(
#        AmphiDEB.M2_complete_ODE_with_loglogistic_TD!; 
#        pmod = p -> begin 
#            p.spc.delta_k_M_mt = 0.33; 
#            return p 
#        end
#        )
#end

@testset "Effect of gamma on shape of the growth trajectory" begin

    maxscale(x) = x ./ maximum(x)

    plt = plot(layout = (2,2), leg = true, size = (800,750), xlabel = ["t" "t" "t_scaled"], ylabel = ["W" "[E_mt]" "W_scaled"])
    plot!(plt, subplot = 4, xaxis = false, yaxis = false, grid = false, xlabel = "", ylabel = "")
    
    p.glb.t_max = 60
    p.glb.pathogen_inoculation_time = Inf
    p.glb.dX_in = [20., 20.]
    p.spc.H_p = 55.

    gamma_values = [0.1, 0.25, 0.5, 0.75, 0.9]

    sims_lrv = DataFrame()

    for gamma in gamma_values
        
        p.spc.gamma = gamma
        sim = AmphiDEB.ODE_simulator(
            p
            );
        
        sim_lrv = @subset(sim, isapprox.(1, :larva, atol = 0.1))
        sim_lrv[!,:t_scaled] = maxscale(sim_lrv.t)
        sim_lrv[!,:W_scaled] = maxscale(sim_lrv.S .+ sim_lrv.E_mt)
        sim_lrv[!,:gamma] .= gamma
        append!(sims_lrv, sim_lrv)

        @df sim plot!(plt, :t, :S .+ :E_mt, label = gamma, subplot = 1)
        @df sim plot!(plt, :t, :E_mt ./ (:S .+ :E_mt), label = gamma, subplot = 2)
        @df sim_lrv plot!(plt, :t_scaled, :W_scaled, label = gamma, subplot = 3)

    end

    # high values of gamma lead to an approximately linear growth curve
    # low values lead to a concave shape of the growth trajectory 
    # if we sort values by gamma for the median time-point, 
    # we should thus get monotonically increasing values for the scaled weight

    sim_midpoint = combine(groupby(sims_lrv, :gamma)) do df

        dist_to_mid = @. abs(df.t_scaled - 0.5)
        df_mid = df[argmin(dist_to_mid),:]

    end

    sort!(sim_midpoint, :W_scaled)
    @test sim_midpoint.gamma == gamma_values

    display(plt)

end

using Distributions

@testset "Randomized parameters" begin
    
    global p = deepcopy(Model1.params)
    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.k_M_emb = truncated(Normal(0.11, 0.011), 0, Inf)
    p.spc.eta_AR = truncated(Normal(0.95, 0.095), 0, 1)
    p.spc.H_p = 50.

    S_max_anl = AmphiDEB.calc_S_max_juv(p.spc)

    @time global sim = @replicates Model1.sim_all(
            p, 
            saveat = 1
            ) 10

    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:dI] = vcat(0, diff(sim.I))
    sim[!,:dI_rel] = sim.dI ./ (sim.S .^(2/3))
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        lineplot(:t, :S), 
        lineplot(:t, :H), 
        lineplot(:t, :E_mt_rel), 
        lineplot(:t, :R),
        lineplot(:t, :adult),
        xrotation = 45
        )
    hline!([p.spc.H_j1], subplot=2, color = :gray, linestyle = :dash)

    display(plt)
    
    # check that all life stage indicators max out close to 1

    @test 0.99 < maximum(sim.embryo) < 1.01
    @test 0.99 < maximum(sim.larva) < 1.01
    @test 0.99 < maximum(sim.metamorph) < 1.01
    @test 0.99 < maximum(sim.juvenile) < 1.01
    @test 0.99 < maximum(sim.adult) < 1.01
    
    # check that the sum of life stage indicators is always approximately 1

    sum_indicators = @. sim.embryo + sim.larva + sim.metamorph + sim.juvenile + sim.adult

    @test unique(isapprox.(1, sum_indicators, atol = 1e-3)) == [true]
end

@testset "gamma = equilbrium reserve  fraction" begin
    p = deepcopy(AmphiDEB.defaultparams)
    p.glb.t_max = 200.
    p.spc.H_j1 = Inf
    gamvals = [0.1, 0.25, 0.5, 0.75, 0.9]
    E_mt_rel_fin = similar(gamvals)

    plt = plot(
        ylabel = ["% reserve" "structure"], 
        xlabel = "t",
        layout = (1,2)
    )

    for (i,gamma) in enumerate(gamvals)
        p.spc.gamma = gamma
        sim = AmphiDEB.ODE_simulator(p)
        sim[!,:E_mt_rel] = sim.E_mt ./ (sim.E_mt .+ sim.S)
        
        @df sim plot!(:t, :E_mt_rel, color = i, subplot = 1)
        hline!([gamma], color = i, linestyle = :dash, subplot = 1)

        @df sim plot!(:t, :S, subplot = 2)

        E_mt_rel_fin[i] = sim.E_mt_rel[end]
    end

    display(plt)

    @test unique(isapprox.(E_mt_rel_fin, gamvals, atol = 0.01)) == [true]
end
