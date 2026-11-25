using Pkg; Pkg.activate("test")
using Aqua
using Revise
using AmphiDEB

Aqua.test_all(AmphiDEB)
