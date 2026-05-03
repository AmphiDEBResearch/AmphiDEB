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

p = ComponentVector(glb=Model1.glb, spc=Model1.spc)
p_ind = Model1.generate_individual_params(p)
u0 = Model1.initialize_statevars(p_ind)
tspan = (0,365)

function du!(du, u, p, t)
    Model1.global!(du.glb, u.glb, p.glb, t)
    Model1.embryo!(du, u, p, t)
end


prob = ODEProblem(du!, u0, tspan, p_ind)
@time sim = solve(prob) |> EcotoxSystems.sol_to_df


using StatsPlots
@df sim plot(:t, :S)

u = sol.u[1]
EcotoxSystems.extract_colnames(u[:ind], :ind)



let u = sol.u[1]
    colnames = []
    for k in keys(u)
        @show k
        push!(colnames, EcotoxSystems.extract_colnames(u[k], k))
    end
    @show colnames
end


using Plots
plot(sol)

sol 

plot(sol.t, sol.yvar)

include("../src/Model1/derivatives.jl")



include("test01_ODE_noeffects.jl")
include("test02_IBM_noeffects.jl")
include("test04_toxicity.jl")
include("test05_pathogens.jl")
include("test07_ODE_temperature.jl")
include("test08_starvation_rules.jl")

include("Aqua.jl")