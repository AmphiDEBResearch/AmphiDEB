occursin("AmphiDEB:jl", pwd()) || cd("./lib/AmphiDEB.jl")
include("boilerplate.jl")
using AmphiDEB.Model1
using AmphiDEB.EcotoxSystems.IBM


@testset "Uninhibited growth" begin
    global p = Model1.defaultparams()

    p.glb.t_max = 450.
    p.glb.dX_in_aq = 1000.
    p.glb.dX_in_ter = 1000.
    p.glb.k_V_aq = 0.
    p.glb.k_V_ter = 0.

    p.glb.N0 = 10

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)

    # these settings are not realistic! only for testing purposes
    p.spc.H_p = 50.  # maturity at puberty
    p.spc.aux.tau_R = 1. # reproduction period
    p.spc.aux.h_S = 0.   # hazard rate caused by loss of structure

    @time sim = IBM.simulate(
        p; 
        global_ode! = Model1.food_dynamics_firstorder!,
        individual_ode! = Model1.individual_ODE!,

        global_rules! = Model1.global_rules!,
        individual_rules! = Model1.individual_rules!,

        init_u_glb = Model1.initialize_global_statevars,
        init_u_ind = Model1.initialize_individual_statevars,
        gen_p_ind = Model1.generate_individual_params,

        N0 = 10, 
        record_individuals = true, 
        saveat = 1/24, 
        dt = 1/24
    )


 
    #@test 2500 <= sim.glb.N[end] <= 3500 # expected abundance after 2 years

    p_glb = @df sim.glb plot(
        :t, :N, 
        xrotation = 45, 
        xlabel = "Time [d]", ylabel = "N", title = "Abundance", 
        leftmargin = 5mm
        )

    p_spc = plot()
    #p_spc = @df sim.spc plot(
    #    groupedlineplot(:t, :S, :cohort, ylabel = "S", title = "Structural mass"), 
    #    groupedlineplot(:t, :cum_repro, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
    #    groupedlineplot(:t, :f_X, :cohort, ylabel = "f(X)", ylim = (0.5, 1.01), title = "Scaled funct. response"), 
    #    groupedlineplot(:t, :H, :cohort, ylabel = "H", title = "Maturity"),
    #    xrotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
    #    bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    #)
    
    plot(p_glb, p_spc, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display

    # check expected order of magnitude of population size after 1.5 years of simulation with defaultparams
    @test 500 < maximum(sim.glb.N) < 1000
end

#p.glb.t_max = 300
#VSCodeServer.@profview_allocs AmphiDEB.IBM_simulator(
#    p; 
#    showinfo = 60,  # update every 30 days
#    saveat = 7, # saving weekly output
#    dt = 1/24 # daily timestep - better to turn down to hourly for proper results
#    )

@testset "Simulation with density-dependence" begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 365. * 3
    p.glb.dX_in = [500., 500.]
    p.glb.k_V = [0.1, 0.1]
    p.glb.N0 = 100

    p.spc.X_emb_int = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(365., 36.5), 0, Inf)
    p.spc.H_p = 50.
    p.spc.K_X_lrv = p.spc.K_X_juv = 20.
    p.spc.h_b = -log(1 - 0.005)

    p.spc.h_S = -log(0.75)
    p.spc.S_rel_crit = 0.5

    @time global sim = AmphiDEB.IBM_simulator(
        p; 
        showinfo = 60, 
        saveat = 1, 
        dt = 1/24, 
        record_individuals = false 
        )

    
    #@test 2500 <= sim.glb.N[end] <= 3500 # expected abundance after 2 years

    p1 = @df sim.glb plot(
        :t, :N, 
        xrotation = 45, 
        xlabel = "Time [d]", ylabel = "N", title = "Total abundance", 
        leftmargin = 5mm, bottommargin = 5mm
        )

    p2 = @df sim.glb plot(
        :t, [:N_emb :N_lrv :N_mt :N_juv :N_ad], 
        label = ["Embryos" "Larvae" "Metamorphs" "Juveniles" "Adults"], 
        xlabel = "Time [d]", ylabel = "N", title = "Abundance per life stage", 
        leg = :outertopright, xrotation = 45, bottommargin = 5mm
    )
    
    plot(p1, p2, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display
end
