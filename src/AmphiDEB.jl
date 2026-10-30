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
The reference model as published in Hansul et al. (2026).
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

"""
An update model version with the following modifications:

- `E_mt` is associated with maintenance costs (assumed identical to structure)
- Investment in `E_mt` starts at fertilization, not birth
- At climax, maturity is re-set to 0
- During climax, ingestion rates are forced to 0
- The rules for mobilizing `E_mt` during climax are the same as the rules for the mobilization of the egg buffer by embryos
"""
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
