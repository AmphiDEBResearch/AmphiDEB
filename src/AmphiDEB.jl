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
    using EcotoxSystems.IBM
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
    include("Model1/rules.jl")
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

end # module AmphiDEB
