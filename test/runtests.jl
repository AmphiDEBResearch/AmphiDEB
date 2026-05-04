
if isinteractive()
    using Pkg; Pkg.activate("lib/AmphiDEB.jl/test")
    using Pkg; Pkg.develop(path="lib/EcotoxSystems.jl")
    using Pkg; Pkg.develop(path="lib/AmphiDEB.jl")
end

using Test, Revise

# TODO: 
#   - update tests to accomodate revised package organization
#   - update IBM configuration to accomodate revised package organization
#       - make sure that switching between life stages remains functional

include("test01_ODE_noeffects.jl")
include("test02_IBM_noeffects.jl")
include("test04_toxicity.jl")
include("test05_pathogens.jl")
include("test07_ODE_temperature.jl")
include("test08_starvation_rules.jl")

include("Aqua.jl")