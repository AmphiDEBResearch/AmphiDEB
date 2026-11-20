using Pkg; Pkg.activate("test")
isinteractive() && Pkg.develop(path="../EcotoxSystems.jl")
Pkg.develop(path=".")

using Chain
using DataFrames 
using DataFramesMeta
using Distributions
using OrdinaryDiffEq
using Plots
using Plots.Measures
using StatsBase
using StatsPlots
using Test

using Revise
using AmphiDEB

default(leg = false, lw = 1.5)
norm(x) = x ./ sum(x)

include("testutils.jl")