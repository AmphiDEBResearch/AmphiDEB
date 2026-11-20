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

- `E_mt` is not implicitly represented before climax and thus simply part of structure
- The model behaviour before climax is therefore identical to the generic DEBkiss model
- This implies that investment in `E_mt` starts at fertilization, not birth and that `E_mt` implies maintenance
- At climax, a fraction `gamma` of structure is converted to `E_mt`
- At climax, maturity is re-set to 0
- During climax
    - Ingestion rates are forced to 0
    - `E_mt` contributes to maintenance and assimilation scaling, since it was derived from structure 
    - Otherwise, the rules for mobilizing `E_mt` during climax are the same as the rules for the mobilization of the egg buffer by embryos
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
