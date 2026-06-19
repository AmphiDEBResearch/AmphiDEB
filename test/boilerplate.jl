using Pkg; Pkg.activate("test")
Pkg.develop(path=".")

using Test
using Plots, StatsPlots, Plots.Measures
using Distributions
using OrdinaryDiffEq
default(leg = false, lw = 1.5)

using Chain
using DataFrames, DataFramesMeta
using StatsBase

using Revise
using AmphiDEB, AmphiDEB.EcotoxSystems