using Pkg; Pkg.activate("test")

using BenchmarkTools

using AmphiDEB.Model2

for _ in 1:10
    p = Model2.defaultparams()
    sim = Model2.sim_all(p; alg = Tsit5())
end

@benchmark Model2.sim_all(p; alg = Tsit5())