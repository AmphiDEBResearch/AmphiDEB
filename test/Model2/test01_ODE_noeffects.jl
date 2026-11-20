using AmphiDEB.Model2
@testset "test default parameters" begin
    alg = Tsit5()

    p = Model2.defaultparams()

    p.glb.t_max = 365 * 2
    p.glb.dX_in_aq = 200. 
    p.glb.dX_in_ter = 200.
    p.glb.food_dynamic = 0.
    
    p.spc.gamma = 0.8
    p.spc.H_p = 55.

    # embryonic stage

    sim_emb, ub, p_ind = Model2.sim_embryo(p; alg = alg)

    @test sim_emb.X_emb[end] ≈ 0 atol = 1e-6
    @test sim_emb.E_mt[end] == 0

    # larval stage

    # simulate until climax
    sim_lrv, uj1 = Model2.sim_larva(p_ind, ub; alg = alg)

    @test sim_lrv.H[end-1] ≈ p.spc.H_j1
    
    H_j1_mem = p_ind.ind.H_j1
    p_ind.ind.H_j1 = Inf

    # simulate longer to check equilibria
    sim_lrv_long, _ = Model2.sim_larva(p_ind, ub; alg = alg, callback = CallbackSet())

    S_max_theo = Model2.get_Smax_lrv(p.spc)
    H_max_theo = Model2.get_Hmax_lrv(p.spc)

    if isinteractive()
        @df sim_lrv_long plot(:t, [:S :H], layout = (1,2))
        hline!([S_max_theo], subplot = 1, color = :gray, linestyle = :dash)
        hline!([H_max_theo], subplot = 2, color = :gray, linestyle = :dash)
    end

    @test minimum(sim_lrv_long.f_X) == 1.
    @test sim_lrv_long.S[end] ≈ S_max_theo atol = 1
    @test sim_lrv_long.H[end] ≈ H_max_theo atol = 1
    @test sim_lrv_long.E_mt[end] == 0

    p_ind.ind.H_j1 = H_j1_mem # reset H_j

    # metamorph stage

    sim_mt, uj2 = Model2.sim_metamorph(
        p_ind, uj1; 
        alg = alg
        )

    @test sim_mt.E_mt[1] ≈ sim_lrv.S[end] * p_ind.ind.gamma atol = 1
    @test sim_mt.E_mt[end] ≈ 0 atol = 1e-6
    
    # juvenile stage

    sim_juv, up = Model2.sim_juvenile(p_ind, uj2; alg = alg)
    @test sim_juv.H[end] ≈ p_ind.ind.H_p

    # adult stage

    sim_ad, u_end = Model2.sim_adult(p_ind, up; alg = alg)
    
    S_max_theo = Model2.get_Smax_ad(p.spc)
    @test minimum(sim_ad.f_X) == 1.
    @test sim_ad.S[end] ≈ S_max_theo atol = 1
    @test sim_ad.E_mt[end] == sim_mt.E_mt[end]
end


#@testset "Randomized parameters" begin
#    
#    global p = deepcopy(Model1.params)
#
#    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
#    p.spc.k_M_emb = truncated(Normal(0.11, 0.011), 0, Inf)
#    p.spc.eta_AR = truncated(Normal(0.95, 0.095), 0, 1)
#    p.spc.H_p = 50.
#
#    S_max_anl = AmphiDEB.calc_S_max_juv(p.spc)
#
#    @time global sim = @replicates Model1.sim_all(
#            p, 
#            saveat = 1
#            ) 10
#
#    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
#    sim[!,:dI] = vcat(0, diff(sim.I))
#    sim[!,:dI_rel] = sim.dI ./ (sim.S .^(2/3))
#    sim[!,:W_tot] = sim.S .+ sim.E_mt 
#    
#    plt = @df sim plot(
#        lineplot(:t, :S), 
#        lineplot(:t, :H), 
#        lineplot(:t, :E_mt_rel), 
#        lineplot(:t, :R),
#        lineplot(:t, :adult),
#        xrotation = 45
#        )
#    hline!([p.spc.H_j1], subplot=2, color = :gray, linestyle = :dash)
#
#    display(plt)
#    
#    # check that all life stage indicators max out close to 1
#
#    @test 0.99 < maximum(sim.embryo) < 1.01
#    @test 0.99 < maximum(sim.larva) < 1.01
#    @test 0.99 < maximum(sim.metamorph) < 1.01
#    @test 0.99 < maximum(sim.juvenile) < 1.01
#    @test 0.99 < maximum(sim.adult) < 1.01
#    
#    # check that the sum of life stage indicators is always approximately 1
#
#    sum_indicators = @. sim.embryo + sim.larva + sim.metamorph + sim.juvenile + sim.adult
#
#    @test unique(isapprox.(1, sum_indicators, atol = 1e-3)) == [true]
#end
