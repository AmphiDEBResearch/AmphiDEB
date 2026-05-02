module AmphiDEB

using PrecompileTools
using EcotoxSystems
import EcotoxSystems: sig, clipneg

using Parameters
using Distributions
using DataStructures
import DataFrames: AbstractDataFrame
using ComponentArrays, StaticArrays
using OrdinaryDiffEq
using StatsBase


include("AmphiDEB_M1/AmphiDEB_M1.jl")
include("AmphiDEB_M1/parameters.jl")
include("AmphiDEB_M1/statevars.jl")
include("AmphiDEB_M1/derivatives.jl")
include("AmphiDEB_M1/rules.jl")
include("AmphiDEB_M1/traits.jl")



include("utils.jl") # various auxiliary functions

# to precompile the model, we simulate the default parameters
@compile_workload begin
    p = deepcopy(params_M1)

    sim = @replicates ODE_simulator(p) 10
    
    p.glb.t_max = 365.
    p.glb.dX_in = [500., 500.]
    p.spc.tau_R = 30.
    
    sim = IBM_simulator(p)
end

end # module AmphiDEB
