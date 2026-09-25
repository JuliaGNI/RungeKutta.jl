using Aqua
using RungeKutta
using Test

Aqua.test_all(RungeKutta; deps_compat = (; broken = true))   # issue #32
