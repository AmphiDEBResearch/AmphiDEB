import ..AmphiDEB.Model2: birth, birth_terminal
import ..AmphiDEB.Model2: puberty, puberty_terminal
import ..AmphiDEB.Model1: aquatic_td_binary_IA_loglogistic

function metamorphosis_condition(u, t, integrator)
    _, _, _, _, y_H, _, _ = aquatic_td_binary_IA_loglogistic(integrator.du, u, integrator.p, t)
    u.ind.H - integrator.p.ind.H_j1 * y_H 
end

function metamorphosis_affect!(integrator)
    integrator.u.ind.embryo = 0.
    integrator.u.ind.larva = 0.
    integrator.u.ind.metamorph = 1.
    integrator.u.ind.juvenile = 0.
    integrator.u.ind.adult = 0.

    integrator.u.ind.H = 0.
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
    save_positions = (true, true)
)

function froglet_emergence_condition(u, t, integrator)
    u.ind.E_mt
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
    nothing;
    save_positions = (true,true)
)

froglet_emergence_terminal = ContinuousCallback(
    froglet_emergence_condition, 
    froglet_emergence_affect_terminal!,
    save_positions = (true,true)
)

"""
Callback condition that hits 0 when structure falls below a minimum value `S_min`.
Callback should be defined to only trigger when `S` downcrosses `S_min`.
We define the condition as `S_min - S`. 
Thus, if the starting `S` is very low (<`S_min`) and `S` crosses `S_min` via regular somatic growth, the callback sees this at downcrossing.
However, if `S` is high and crosses `S_min` via shrinking, the callback sees this at upcrossing and should be triggered (e.g. to terminate the integration). 
"""
function min_structure_condition(u, t, integrator)
    p.ind.aux.S_min - u.ind.S
end
