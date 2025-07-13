"""
    AmphiDEB_global_rules!(m)

Global rule-based portion of the AmphiDEB IBM. 
"""
function AmphiDEB_global_rules!(m)
    m.u.glb.N = length(m.individuals) # tracking population size
    m.u.glb.X = max.(0, m.u.glb.X) # HOTFIX : negative resource abundances can cause chaos
    
    # record counts per life stage
    m.u.glb.N_emb = sum(map(x -> isapprox(1, x.u.ind.embryo, atol = 0.01), m.individuals))
    m.u.glb.N_lrv = sum(map(x -> isapprox(1, x.u.ind.larva, atol = 0.01), m.individuals))
    m.u.glb.N_mt = sum(map(x -> isapprox(1, x.u.ind.metamorph, atol = 0.01), m.individuals))
    m.u.glb.N_juv = sum(map(x -> isapprox(1, x.u.ind.juvenile, atol = 0.01), m.individuals))
    m.u.glb.N_ad = sum(map(x -> isapprox(1, x.u.ind.adult, atol = 0.01), m.individuals))
    
    m.u.glb.N_emb = sum(isapprox(1, x.u.ind.embryo; atol=0.01) for x in m.individuals)
    m.u.glb.N_lrv = sum(isapprox(1, x.u.ind.larva;  atol=0.01) for x in m.individuals)
    m.u.glb.N_mt  = sum(isapprox(1, x.u.ind.metamorph; atol=0.01) for x in m.individuals)
    m.u.glb.N_juv = sum(isapprox(1, x.u.ind.juvenile; atol=0.01) for x in m.individuals)
    m.u.glb.N_ad  = sum(isapprox(1, x.u.ind.adult;    atol=0.01) for x in m.individuals)
end
