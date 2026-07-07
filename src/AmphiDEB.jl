module AmphiDEB

using PrecompileTools
using EcotoxSystems

using Parameters
using Distributions
using DataStructures
import DataFrames: AbstractDataFrame
using ComponentArrays, StaticArrays # TODO: static arrays can probably go
using OrdinaryDiffEq
using StatsBase

include("utils.jl") # various auxiliary functions


module FittingUtils

    using ComponentArrays
    using DataFrames
    using Distributions
    using Chain

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
    include("Model1/simulation.jl")
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

module Model2b
    using EcotoxSystems, EcotoxModelFitting
    import ..FittingUtils

    using Parameters
    using Distributions
    using DataStructures
    import DataFrames: AbstractDataFrame
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase

    include("Model2b/parameters.jl")
    include("Model2b/statevars.jl")
    include("Model2b/callbacks.jl")
    include("Model2b/derivatives.jl")

end

module Model2c
    using EcotoxSystems, EcotoxModelFitting
    import ..FittingUtils

    using Parameters
    using Distributions
    using DataStructures
    import DataFrames: AbstractDataFrame
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase

    include("Model2c/parameters.jl")
    include("Model2c/statevars.jl")
    include("Model2c/callbacks.jl")
    include("Model2c/derivatives.jl")
end


module Model2d
    using EcotoxSystems, EcotoxModelFitting
    import ..FittingUtils

    using Parameters
    using Distributions
    using DataStructures
    import DataFrames: AbstractDataFrame
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase
    using DataFrames

    include("Model2d/parameters.jl")
    include("Model2d/statevars.jl")
    include("Model2d/callbacks.jl")
    include("Model2d/derivatives.jl")
    include("Model2d/traits.jl")
end

module Model3
    using EcotoxSystems, EcotoxModelFitting
    import ..FittingUtils
    import ..AmphiDEB: sabs, smax

    using Parameters
    using Distributions
    using DataStructures
    import DataFrames: AbstractDataFrame
    using ComponentArrays, StaticArrays
    using OrdinaryDiffEq
    using StatsBase
    using DataFrames

    include("Model3/parameters.jl")
    include("Model3/statevars.jl")
    include("Model3/callbacks.jl")
    include("Model3/derivatives.jl")
    include("Model3/traits.jl")
end


end # module AmphiDEB
