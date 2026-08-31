# Runge-Kutta Methods

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/RungeKutta.jl/stable)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliagni.github.io/RungeKutta.jl/latest)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/R/RungeKutta.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/R/RungeKutta.html)
[![Build Status](https://github.com/JuliaGNI/RungeKutta.jl/workflows/CI/badge.svg)](https://github.com/JuliaGNI/RungeKutta.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaGNI/RungeKutta.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/RungeKutta.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.4294923.svg)](https://doi.org/10.5281/zenodo.4294923)

This package collects Runge-Kutta tableaus and provides diagnostics to analyze them.
It implements algorithms for the computation of Gauss, Lobatto and Radau tableaus with arbitrary numbers of stages as well as tabulated coefficients for various explicit, diagonally implicit and fully implicit methods.
All tableaus can be retrieved in arbitrary precision or symbolically.

## Installation

*RungeKutta.jl* and all of its dependencies can be installed via the Julia REPL by typing 
```julia
]add RungeKutta
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

- *explicit*: ExplicitEuler/ForwardEuler, ExplicitMidpoint, Heun2, Heun3, Ralston2, Ralston3, Runge2, Kutta3, RK21, RK22, RK31, RK32, RK4, RK41, RK42, RK416, RK438, RK5, SSPRK2, SSPRK3

- *diagonally implicit*: CrankNicolson, Crouzeix, KraaijevangerSpijker, QinZhang

- *fully implicit*: ImplicitEuler/BackwardEuler, ImplicitMidpoint, IRK3, SRK3

In addition there exist functions to compute Gauss, Lobatto and Radau tableaus with an arbitrary number of stages s:

- `TableauGauss(s)`

- `TableauLobattoIIIA(s)`, `TableauLobattoIIIB(s)`, `TableauLobattoIIIC(s)`, `TableauLobattoIIIC̄(s)`, `TableauLobattoIIID(s)`, `TableauLobattoIIIE(s)`, `TableauLobattoIIIF(s)`, `TableauLobattoIIIG(s)`

- `TableauRadauIA(s)`, `TableauRadauIIA(s)`

All constructors take an optional type argument `T`, as in `TableauExplicitMidpoint(T)` or `TableauGauss(T,s)`. The default type is `Float64`, but it can be set to other number types if needed, e.g., to `Float32` for single precision or to the `Dec128` type from [DecFP.jl](https://github.com/JuliaMath/DecFP.jl) for quadruple precision.
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
* `satisfies_simplifying_assumption_b(tab, σ=tab.s)`
* `satisfies_simplifying_assumption_c(tab, σ=tab.s)`
* `satisfies_simplifying_assumption_d(tab, σ=tab.s)`

This list is expected to grow in the near future.


## References

If you use RungeKutta.jl in your work, please consider citing it by

```
@misc{Kraus:2020:RungeKutta,
  title={RungeKutta.jl: Runge-Kutta Methods in Julia},
  author={Kraus, Michael},
  year={2020},
  howpublished={\url{https://github.com/JuliaGNI/RungeKutta.jl}},
  doi={10.5281/zenodo.4294923}
}
```


## Development

### Git hooks

Two hooks live in `.githooks`. They are **not active in a fresh clone** — `core.hooksPath` is local
configuration and does not travel with a push — so enable them once per clone:

```sh
git config core.hooksPath .githooks
```

**`pre-commit`** acts on **staged `.jl` files only**, and exits immediately when a commit stages
none, so a documentation- or workflow-only commit is not slowed down by it:

- **JuliaFormatter `--check`**, honouring this repository's own `.JuliaFormatter.toml` — **blocks**
  the commit. Formatting is mechanical and always fixable.
- **`fatou lint`**, when `fatou` is installed — **advisory only**, and deliberately so: its
  `unused-import` rule does not follow `include`, so it flags the load-bearing imports of every
  module file.
- **`using <Package>`**, which catches a syntax error or a broken `include` — **blocks**.

**`pre-push`** runs the full test suite with `--check-bounds=auto`, but **only when pushing to
`main` or `master`**; a topic branch is left to CI. It prints nothing for **10–30 minutes**, which
looks exactly like a network hang and is not one. If you do interrupt it, check for an orphaned
Julia process that the killed hook left behind.

Either hook can be bypassed for a single command with `--no-verify`, for a change you know it does
not apply to:

```sh
git commit --no-verify
git push --no-verify
```

The hooks are generated from one shared copy and are byte-identical across the related
repositories, so edit them there rather than here — a local edit is silently undone by the next
install.
