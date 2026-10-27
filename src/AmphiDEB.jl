module AmphiDEB

using PrecompileTools
using EcotoxSystems

using Parameters
using Distributions
using DataStructures
using DataFrames
using ComponentArrays, StaticArrays # TODO: static arrays can probably go
using OrdinaryDiffEq
using StatsBase

include("utils.jl") # various auxiliary functions

"""
The reference model.
"""
module Model1

    using EcotoxSystems, EcotoxModelFitting
    using EcotoxSystems.IBM

    using Parameters
    using Distributions
    using DataFrames
    using DataStructures
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase

    include("Model1/parameters.jl")
    include("Model1/statevars.jl")
    include("Model1/rules.jl")
    include("Model1/callbacks.jl")
    include("Model1/derivatives.jl")
    include("Model1/simulation.jl")
    include("Model1/traits.jl")
end

# TODO: add a separate `f` for metamorphs to correct for differences in ingestion rate 
#       vs. tissue mobilization rate 
module Model2

using EcotoxSystems, EcotoxModelFitting
    using EcotoxSystems.IBM

    using Parameters
    using Distributions
    using DataStructures
    import DataFrames: AbstractDataFrame, DataFrameRow
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase

    include("Model2/parameters.jl")
    include("Model2/statevars.jl")
    include("Model2/rules.jl")
    include("Model2/callbacks.jl")
    include("Model2/derivatives.jl")
    include("Model2/simulation.jl")
    include("Model2/traits.jl")

end


end # module AmphiDEB
