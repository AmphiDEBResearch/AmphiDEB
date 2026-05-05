module AmphiDEB

using PrecompileTools
using EcotoxSystems

using Parameters
using Distributions
using DataStructures
import DataFrames: AbstractDataFrame
using ComponentArrays, StaticArrays
using OrdinaryDiffEq
using StatsBase


module FittingUtils
    include("FittingUtils/losses.jl")
    include("FittingUtils/utils.jl")
    include("FittingUtils/traits.jl")
end



module Model1

    using EcotoxSystems

    using Parameters
    using Distributions
    using DataStructures
    import DataFrames: AbstractDataFrame
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase
    include("Model1/parameters.jl")
    include("Model1/statevars.jl")
    include("Model1/callbacks.jl")
    include("Model1/derivatives.jl")
    include("Model1/fit_larvae_fulldata.jl")
end


# for backwards-compatability:

import .Model1 
defaultparams = Model1.params # old "defaultparams" is an alias for "Model1.params"
const ODE_simulator = Model1.sim_all # old "ODE_simulator" is an alias for the now "Model1.sim_all"

include("individual_rules.jl") # individual rule-based component
include("global_rules.jl") # global rule-based component
#include("simulators.jl") # functions to run run simulations
include("traits.jl") 

include("utils.jl") # various auxiliary functions

# to precompile the model, we simulate the default parameters
@compile_workload begin
    p = deepcopy(defaultparams)

    sim = @replicates ODE_simulator(p) 10
end



end # module AmphiDEB
