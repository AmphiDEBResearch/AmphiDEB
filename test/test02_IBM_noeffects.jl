occursin("AmphiDEB.jl", pwd()) || cd("./lib/AmphiDEB.jl")
include("boilerplate.jl")
using AmphiDEB.Model1
using AmphiDEB.EcotoxSystems.IBM

Model1.defaultparams()

# FIXME: why is N_emb, N_lrv...not updated correctly?

begin
    global p = Model1.defaultparams()


    # these settings are not realistic - only for testing purposes

    p.glb.t_max = 10 * 365
    p.glb.dX_in_aq = 1000.
    p.glb.dX_in_ter = 10_000.
    p.glb.k_V_aq = 0.1
    p.glb.k_V_ter = 0.1

    p.glb.N0 = 10

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)

    p.spc.H_p = 50.  # maturity at puberty
    p.spc.aux.tau_R = Truncated(Normal(365., 36.5), 0, Inf) # reproduction period
    p.spc.aux.h_S = 0.1 # hazard rate caused by loss of structure
    p.spc.aux.S_rel_crit = 0.5
    p.spc.K_X_lrv = 100.
    p.spc.K_X_juv = 100.

    @time sim = IBM.simulate(
        p; 
        global_ode! = Model1.food_dynamics_firstorder!,
        individual_ode! = Model1.individual_ODE!,

        global_rules! = Model1.global_rules!,
        individual_rules! = Model1.individual_rules!,

        init_u_glb = Model1.initialize_global_statevars,
        init_u_ind = Model1.initialize_individual_statevars,
        gen_p_ind = Model1.generate_individual_params,

        N0 = 100,
        record_individuals = false, 
        saveat = 1, 
        dt = 1/24
    )

    p_glb = @df sim.glb plot(
        plot(
            :t, [:N :N_emb :N_lrv :N_mt :N_juv :N_ad], 
            leg = :outertopright, 
            label = ["total" "embryos" "larvae" "metamorphs" "juveniles" "adults"],
            legendtitle = "life stage",
            leftmargin = 5mm, 
            ylabel = "N [#]", title = "Abundance", 
            ), 
        plot(
            :t, [:X_aq :X_ter], 
            leg = :outertopright, label = ["aquatic" "terrestric"], 
            legendtitle = "habitat",
            ylabel = "X [mg]", title = "Food"
            ), 
        size = (1000,450), thickness_scaling = 1.2, 
        bottommargin = 5mm,
        xlabel = "Time [d]", 
    )

    p_glb |> display
    #@transform! sim.spc :fX = map(x->x.fX, :aux)
    #@df @subset(sim.spc, :juvenile .== 1) lineplot(:t, :fX)
end


#@testset "Example simulation with IBM" begin
#   # TODO: implement proper test
#end



