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


#module FittingUtils
#
#    using ComponentArrays
#    using DataFrames
#    using Distributions
#    using Chain
#
#    include("FittingUtils/losses.jl")
#    include("FittingUtils/utils.jl")
#    include("FittingUtils/traits.jl")
#end

"""
The reference model.
"""
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
    include("Model1/traits.jl")
end

module Model2

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

    include("Model2/parameters.jl")
    include("Model2/statevars.jl")
    include("Model2/rules.jl")
    include("Model2/callbacks.jl")
    include("Model2/derivatives.jl")
    include("Model2/simulation.jl")

end

#"""
#`Model1b` is identical to `Model1`, except that maturity is re-set at climax (`j1`, GS 42).
#"""
#module Model1b 
#    using EcotoxSystems, EcotoxModelFitting
#    using EcotoxSystems.IBM
#    import ..FittingUtils
#
#    using Parameters
#    using Distributions
#    using DataStructures
#    import DataFrames: AbstractDataFrame
#    using ComponentArrays, StaticArrays
#    using OrdinaryDiffEq
#    using StatsBase
#
#    include("Model1b/parameters.jl")
#    include("Model1b/statevars.jl")
#    include("Model1b/rules.jl")
#    include("Model1b/callbacks.jl")
#    include("Model1b/derivatives.jl")
#    include("Model1b/simulation.jl")
#
#end
#
#"""
#`Model2` forces ingestion rates to 0 during climax. 
#Mobilization of the buffer follows a first-order rate constant.
#Froglet emergence occurs when the buffer is used up, up to a rest given by the auxiliary parameter `ϵ_Emt`.
#"""
#module Model2
#    using EcotoxSystems, EcotoxModelFitting
#    import ..FittingUtils
#
#    using Parameters
#    using Distributions
#    using DataStructures
#    import DataFrames: AbstractDataFrame
#    using ComponentArrays, StaticArrays
#    using OrdinaryDiffEq
#    using StatsBase
#
#    include("Model2/parameters.jl")
#    include("Model2/statevars.jl")
#    include("Model2/callbacks.jl")
#    include("Model2/derivatives.jl")
#
#end
#
#"""
#`Model2b` is identical to `Model2`, except that maturity is re-set at climax.
#"""
#module Model2b
#    using EcotoxSystems, EcotoxModelFitting
#    import ..FittingUtils
#
#    using Parameters
#    using Distributions
#    using DataStructures
#    import DataFrames: AbstractDataFrame
#    using ComponentArrays, StaticArrays
#    using OrdinaryDiffEq
#    using StatsBase
#
#    include("Model2b/parameters.jl")
#    include("Model2b/statevars.jl")
#    include("Model2b/callbacks.jl")
#    include("Model2b/derivatives.jl")
#
#end
#
#"""
#`Model2c` includes maintenance costs associated with `E_mt`.
#"""
#module Model2c
#    using EcotoxSystems, EcotoxModelFitting
#    import ..FittingUtils
#
#    using Parameters
#    using Distributions
#    using DataStructures
#    import DataFrames: AbstractDataFrame
#    using ComponentArrays, StaticArrays
#    using OrdinaryDiffEq
#    using StatsBase
#
#    include("Model2c/parameters.jl")
#    include("Model2c/statevars.jl")
#    include("Model2c/callbacks.jl")
#    include("Model2c/derivatives.jl")
#end
#
#"""
#`Model2d`
#"""
#module Model2d
#    using EcotoxSystems, EcotoxModelFitting
#    import ..FittingUtils
#
#    using Parameters
#    using Distributions
#    using DataStructures
#    import DataFrames: AbstractDataFrame
#    using ComponentArrays, StaticArrays
#    using OrdinaryDiffEq
#    using StatsBase
#    using DataFrames
#
#    include("Model2d/parameters.jl")
#    include("Model2d/statevars.jl")
#    include("Model2d/callbacks.jl")
#    include("Model2d/derivatives.jl")
#    include("Model2d/traits.jl")
#end
#
#"""
#In `Model3`, froglet emergence is trigerred by a maturity threshold `H_j1`.
#At emergence, froglets convert remaining reserve buffer to structure, so that froglets can continue to function like a standard DEBkiss animal. 
#
#Investment in `E_mt` starts at fertilization, not birth. 
#`E_mt` is associated with specific maintenance costs.
#"""
#module Model3
#    using EcotoxSystems, EcotoxModelFitting
#    import ..FittingUtils
#    import ..AmphiDEB: sabs, smax
#
#    using Parameters
#    using Distributions
#    using DataStructures
#    import DataFrames: AbstractDataFrame
#    using ComponentArrays, StaticArrays
#    using OrdinaryDiffEq
#    using StatsBase
#    using DataFrames
#
#    include("Model3/parameters.jl")
#    include("Model3/statevars.jl")
#    include("Model3/callbacks.jl")
#    include("Model3/derivatives.jl")
#    include("Model3/traits.jl")
#end
#
#
#"""
#Model 4 further simplifies the model structure by removing the state variable E_mt altogether.
#
#Metamorphs mimick the embryonic phase more strictly through the following assumptions:
#
#- Maturity is re-set at climax
#- At climax, a fraction γW_S of structural weight is converted to buffer
#- The buffer has the same function as the egg buffer for embryos. 
#- Following the default starvation rules, metamorphs may use the buffer to cover maintenance.
#"""
#module Model4
#
#
#end

end # module AmphiDEB
