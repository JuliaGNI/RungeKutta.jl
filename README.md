# Runge Kutta

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/RungeKutta.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliagni.github.io/RungeKutta.jl/dev)
[![Build Status](https://github.com/JuliaGNI/RungeKutta.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/RungeKutta.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/RungeKutta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/RungeKutta.jl)

This package collects Runge-Kutta tableaus and provides diagnostics to analyze them.
It implements algorithms for the computation of Gauss, Lobatto and Radau tableaus with arbitrary numbers of stages as well as tabulated coefficients for various explicit, diagonally implicit and fully implicit methods.
All tableaus can be retrieved in arbitrary precision.

## Installation

*RungeKutta.jl* and all of its dependencies can be installed via the Julia REPL by typing 
```julia
]add https://github.com/JuliaGNI/RungeKutta.jl
```

## Basic Usage

After loading the Runge-Kutta module by
```julia
julia> using RungeKutta
```
a `Tableau` can be created by calling any one of the provided constructors, for example
```julia
julia> TableauExplicitMidpoint()

Runge-Kutta Tableau explicit_midpoint with 2 stages and order 2:

 0.0 │ 0.0  0.0
 0.5 │ 0.5  0.0
─────┼──────────
     │ 0.0  1.0
```

The `Tableau` type has the following fields
- `name` is a descriptive name of the tableau,
- `o` the order of the method,
- `s` the number of stages,
- `a` the coefficients,
- `b` the weights,
- `c` the nodes.

The following tableaus are implemented (prepend `Tableau` to the name to call the respective constructor):

- *explicit*: ExplicitEuler/ForwardEuler, ExplicitMidpoint, Heun2, Heun3, Ralston2, Ralston3, Runge/Runge2, Kutta/Kutta3, RK4/RK416, RK438, SSPRK3

- *diagonally implicit*: KraaijevangerSpijker, QinZhang, Crouzeix

- *fully implicit*: ImplicitEuler/BackwardEuler, ImplicitMidpoint, CrankNicolson, SRK3

In addition there exist functions to compute Gauss, Lobatto and Radau tableaus with an arbitrary number of stages s:

- `TableauGauss(s)`

- `TableauLobattoIIIA(s)`, `TableauLobattoIIIB(s)`, `TableauLobattoIIIC(s)`, `TableauLobattoIIIC̄(s)`, `TableauLobattoIIID(s)`, `TableauLobattoIIIE(s)`, `TableauLobattoIIIF(s)`, `TableauLobattoIIIG(s)`

- `TableauRadauIIA(s)`

All constructors take an optional type argument `T`, as in `TableauExplicitMidpoint(T)` or `TableauGauss(s,T)`. The default type is `Float64`, but it can be set to other number types if needed, e.g., to `Float32` for single precision or to the `Dec128` type from [DecFP.jl](https://github.com/JuliaMath/DecFP.jl) for quadruple precision.
Internally, all tableaus are computed using `BigFloat`, providing high-accuracy coefficients as they are required for simulations in quadruple or higher precision. The internal precision can be set via `setprecision(40)`, cf. the [Julia Manual](https://docs.julialang.org/en/v1/) on [Arbitrary Precision Arithmetic](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Arbitrary-Precision-Arithmetic).

## Custom Tableaus

If required, it is straight-forward to create a custom tableau.
The tableau of Heun's method, for example, is defined as follows:
```julia
a = [[0.0  0.0]
     [1.0  0.0]]
b =  [0.5, 0.5]
c =  [0.0, 1.0]
o = 2

tab = Tableau(:heun, o, a, b, c)
```

## Diagnostics

Currently, diagnostic functions for checking symmetry, symplecticity and the so-called simplifying assumptions are implemented:

* `issymmetric(tab)`
* `issymplectic(tab)`
* `satisfies_simplifying_assumption_b(tab, σ=s)`
* `satisfies_simplifying_assumption_c(tab, σ=s)`
* `satisfies_simplifying_assumption_d(tab, σ=s)`

This list is expected to grow in the near future.
