import ..AmphiDEB.Model2: birth, birth_terminal
import ..AmphiDEB.Model2: puberty, puberty_terminal
import ..AmphiDEB.Model2: froglet_emergence, froglet_emergence_terminal
import ..AmphiDEB.Model1: aquatic_td_binary_IA_loglogistic

# everything identical to Model2 except for metamorphosis callback

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