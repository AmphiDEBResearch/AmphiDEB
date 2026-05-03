if !isinteractive()
    using Pkg; Pkg.activate("test")
end

if isinteractive()
    using Pkg; Pkg.activate("lib/AmphiDEB.jl/test")
end

using Pkg; Pkg.develop(path="lib/AmphiDEB.jl")

Pkg.instantiate()
Pkg.resolve()

using Test
using Plots, StatsPlots, Plots.Measures
using Distributions
using OrdinaryDiffEq
default(leg = false, lw = 1.5)

using Chain
using DataFrames, DataFramesMeta
using StatsBase

using Revise

using AmphiDEB
using EcotoxSystems

import EcotoxSystems: sig
import EcotoxSystems: constrmvec

include("test01_ODE_noeffects.jl")
include("test02_IBM_noeffects.jl")
include("test04_toxicity.jl")
include("test05_pathogens.jl")
include("test07_ODE_temperature.jl")
include("test08_starvation_rules.jl")

include("Aqua.jl")