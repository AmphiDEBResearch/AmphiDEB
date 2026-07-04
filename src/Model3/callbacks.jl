birth_condition(u, t, integrator) = u.ind.X_emb
function birth_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 1.
    integrator.u.ind.metamorph = 0.
    integrator.u.ind.juvenile = 0.
    integrator.u.ind.adult = 0.
end

function birth_affect_terminal!(integrator)
    birth_affect!(integrator)
    terminate!(integrator)
end

birth = ContinuousCallback(
    birth_condition, 
    birth_affect!, 
    nothing, # we need `neg_affect! = nothing` to tell the solver that the affect should only occur for upcrossings (condition function switches from negative to positive) 
    save_positions = (true,true)
    )

birth_terminal = ContinuousCallback(
    birth_condition,
    birth_affect_terminal!,
    save_positions = (true,true)
)

function metamorphosis_condition(u, t, integrator)
    _, _, _, _, y_H, _, _ = aquatic_td_binary_IA_loglogistic(integrator.du, u, integrator.p, t)
    u.ind.E_mt - integrator.p.ind.H_j1 * y_H 
end

function metamorphosis_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 0.
    integrator.u.ind.metamorph = 1.
    integrator.u.ind.juvenile = 0.
    integrator.u.ind.adult = 0.
end

function metamorphosis_affect_terminal!(integrator)
    metamorphosis_affect!(integrator)
    terminate!(integrator)
end

metamorphosis = ContinuousCallback(
    metamorphosis_condition, 
    metamorphosis_affect!,
    save_positions = (true,true)
)

metamorphosis_terminal = ContinuousCallback(
    metamorphosis_condition, 
    metamorphosis_affect_terminal!,
    save_positions = (true,true)
)

function froglet_emergence_condition(u, t, integrator)
    integrator.p.ind.H_j2 - u.ind.E_mt
end

function froglet_emergence_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 0.
    integrator.u.ind.metamorph = 0.
    integrator.u.ind.juvenile = 1.
    integrator.u.ind.adult = 0.

    # at emergence, froglets convert remaining reserve buffer to structure
    integrator.u.ind.S += integrator.u.ind.E_mt
    integrator.u.ind.E_mt = 0.
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

puberty_condition(u, t, integrator) = u.ind.H - integrator.p.ind.H_p

function puberty_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 0.
    integrator.u.ind.metamorph = 0.
    integrator.u.ind.juvenile = 0.
    integrator.u.ind.adult = 1.
end

function puberty_affect_terminal!(integrator)
    puberty_affect!(integrator)
    terminate!(integrator)
end

puberty = ContinuousCallback(
    puberty_condition, 
    puberty_affect!, 
    nothing,
    save_positions = (true,true)
    )
puberty_terminal = ContinuousCallback(
    puberty_condition, 
    puberty_affect_terminal!, 
    nothing,
    save_positions = (true,true)
    )

