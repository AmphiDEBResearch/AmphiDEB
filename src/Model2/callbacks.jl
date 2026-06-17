import ..AmphiDEB.Model1: birth, birth_terminal, metamorphosis, metamorphosis_terminal, puberty, puberty_terminal

function froglet_emergence_condition(u, t, integrator)
    u.ind.E_mt - integrator.p.ind.aux.ϵ_Emt
end

function froglet_emergence_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 0.
    integrator.u.ind.metamorph = 0.
    integrator.u.ind.juvenile = 1.
    integrator.u.ind.adult = 0.
end

function froglet_emergence_affect_terminal!(integrator)
    froglet_emergence_affect!(integrator)
    terminate!(integrator)
end

froglet_emergence = ContinuousCallback(
    froglet_emergence_condition, 
    froglet_emergence_affect!,
    save_positions = (true,true)
)

froglet_emergence_terminal = ContinuousCallback(
    froglet_emergence_condition, 
    froglet_emergence_affect_terminal!,
    save_positions = (true,true)
)
