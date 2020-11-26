```@meta
CurrentModule = RungeKutta
```

# Runge Kutta

[![Build Status](https://github.com/JuliaGNI/RungeKutta.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/RungeKutta.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/RungeKutta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/RungeKutta.jl)

This package collects Runge-Kutta tableaus and diagnostics to analyze them.
All tableaus are either computed or evaluated using `BigFloat` leading to high-precision coefficients. This is important for simulations in quadruple or higher precision.
We provide algorithms for the computation of Gauss, Lobatto and Radau tableaus with arbitrary numbers of stages as well as tabulated tableaus for various explicit, diagonally and fully implicit methods.

