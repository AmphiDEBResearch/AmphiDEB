if !isinteractive()
    using Pkg; Pkg.activate("test")
end

if isinteractive()
    using Pkg; Pkg.activate("lib/AmphiDEB.jl/test")
    using Pkg; Pkg.develop(path="lib/EcotoxSystems.jl")
    using Pkg; Pkg.develop(path="lib/AmphiDEB.jl")
end


using Revise
using AmphiDEB.Model1

using AmphiDEB.Parameters, AmphiDEB.ComponentArrays, AmphiDEB.Distributions
using EcotoxSystems
using OrdinaryDiffEq

using StatsPlots


sim_emb, u0lrv, p_ind = Model1.sim_embryo(p)

u0lrv

begin
    p = ComponentVector(glb=Model1.glb, spc=Model1.spc)
    p_ind = Model1.generate_individual_params(p)
    u0 = Model1.initialize_statevars(p_ind)
    tspan = (0,365)


    prob = ODEProblem(Model1.du!, u0, tspan, p_ind)
    @time sim = solve(prob, callback = Model1.callback_set) |> EcotoxSystems.sol_to_df

    @df sim plot(
        plot(:t, :S, ylabel = "S"), 
        plot(:t, :H, ylabel = "H"), 
        plot(:t, :X_emb, ylabel = "X_emb"),
        plot(:t, [:is_embryo :is_larva :is_metamorph])
    )

    hline!(p.spc.H_j1, subplot = 2)
end


@df sim plot(
    plot(:t; :is_embryo),
    plot(:t, :is_larva), 
    plot(:t, :is_metamorph)
)

sim.is_embryo |> unique
sim.is_larva |> unique
sim.is_metamorph |> unique


include("../src/Model1/derivatives.jl")



include("test01_ODE_noeffects.jl")
include("test02_IBM_noeffects.jl")
include("test04_toxicity.jl")
include("test05_pathogens.jl")
include("test07_ODE_temperature.jl")
include("test08_starvation_rules.jl")

include("Aqua.jl")