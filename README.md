# Runge Kutta

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/RungeKutta.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliagni.github.io/RungeKutta.jl/dev)
[![Build Status](https://github.com/JuliaGNI/RungeKutta.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/RungeKutta.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/RungeKutta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/RungeKutta.jl)

This package collects Runge-Kutta tableaus and provides diagnostics to analyze them.
Most of the functionality is extracted from [GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl).
All tableaus are either computed or evaluated using `BigFloat` leading to high-precision coefficients. This is important for simulations in quadruple or higher precision.
We provide algorithms for the computation of Gauss, Lobatto and Radau tableaus with arbitrary numbers of stages as well as tabulated tableaus for various explicit, diagonally and fully implicit methods.

## Installation

*RungeKutta.jl* and all of its dependencies can be installed via the Julia REPL by typing 
```
]add RungeKutta
```
