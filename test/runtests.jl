cd("./lib/AmphiDEB.jl")
include("boilerplate.jl")

# TODO: Tests for Model 1 probably need to be updated...

#include("Model1/test01_ODE_noeffects.jl")
#include("Model1/test02_IBM_noeffects.jl")
#include("Model1/test04_toxicity.jl")
#include("Model1/test05_pathogens.jl")
#include("Model1/test07_ODE_temperature.jl")
#include("Model1/test08_starvation_rules.jl")

include("Model2/test01_ODE_noeffects.jl")
# TODO: remaining tests for Model 2 need to be updated...

include("Aqua.jl")