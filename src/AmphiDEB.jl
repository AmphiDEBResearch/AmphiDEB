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

    include("FittingUtils/losses.jl")
    include("FittingUtils/utils.jl")
    include("FittingUtils/traits.jl")
end

module Model1

    using EcotoxSystems, EcotoxModelFitting
    import ..FittingUtils

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
end

module Model2
    using EcotoxSystems, EcotoxModelFitting
    import ..FittingUtils

    using Parameters
    using Distributions
    using DataStructures
    import DataFrames: AbstractDataFrame
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase

    include("Model2/parameters.jl")
    include("Model2/statevars.jl")
    include("Model2/callbacks.jl")
    include("Model2/derivatives.jl")

end


## for backwards-compatability:
#
#import .Model1 
#defaultparams = Model1.params # old "defaultparams" is an alias for "Model1.params"
#const ODE_simulator = Model1.sim_all # old "ODE_simulator" is an alias for the now "Model1.sim_all"
#
#include("individual_rules.jl") # individual rule-based component
#include("global_rules.jl") # global rule-based component
##include("simulators.jl") # functions to run run simulations
#include("traits.jl") 
#
#include("utils.jl") # various auxiliary functions

# to precompile the model, we simulate the default parameters
#@compile_workload begin
#    sim = replicates(p->ODE_simulator(p), params, 10)
#end

end # module AmphiDEB
