# Tableau

The [`Tableau`](@ref) type is the main data structure in RungeKutta.jl.
It holds all the coeffieints and information on a Runge-Kutta method in the form of a so-called Butcher tableau

```math
\begin{array}{c|c}
c & a     \\
\hline
  & b^{T} \\
\end{array}
=
\begin{array}{c|cccc}
c_{1}  & a_{11} & a_{12} & \dots & a_{1s} \\
c_{2}  & a_{21} & a_{22} & \dots & a_{2s} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
c_{s}  & a_{s1} & a_{s2} & \dots & a_{ss} \\
\hline
       & b_{1}  & b_{2}  & \dots & b_{s}  \\
\end{array}
```

The [`Tableau`](@ref) type has the following fields:
- `name` is a descriptive name of the tableau,
- `o` the order of the method,
- `s` the number of stages,
- `a` the coefficients,
- `b` the weights,
- `c` the nodes.

Although one of the main goals of RungeKutta.jl is to provide as many known tableaus as possible (see below), it is straight-forward to create a custom tableau.
The tableau of Heun's method, for example, is defined as follows:
```julia
a = [[0.0  0.0]
     [1.0  0.0]]
b =  [0.5, 0.5]
c =  [0.0, 1.0]
o = 2

tab = Tableau(:heun, o, a, b, c)
```

There exist several constructors
```julia
Tableau{T}(name, o, s, a, b, c)
Tableau{T}(name, o, a, b, c)
Tableau(name,o, s, a, b, c)
Tableau(name,o, a, b, c)
```
The meaning of the arguments corresponds to the fields as explained above and `T` is the datatype of the coefficient arrays.
If `s` is ommitted it is inferred from the length of `c`. All coefficient arrays are checked for compatible sizes.
If `T` is not specified, the datatype is inferred from the type of `a`, `b` and `c`, which are then assumed to be identical.
If `T` is specified and does not correspond to the types of  `a`, `b` and `c`, the arrays are converted accordingly.


## Tabulated Coefficients

RungeKutta.jl provides tabulated coefficients for various explicit and implicit (both diagonally and fully implicit) Runge-Kutta methods, in particular

| Function and Aliases                                           | Stages | Order |
|:---------------------------------------------------------------|:-------|:------|
| **Explicit Methods**                                           |        |       |
| [`TableauExplicitEuler`](@ref), [`TableauForwardEuler`](@ref)  | 1      | 1     |
| [`TableauExplicitMidpoint`](@ref)                              | 2      | 2     |
| [`TableauHeun2`](@ref)                                         | 2      | 2     |
| [`TableauHeun3`](@ref)                                         | 3      | 3     |
| [`TableauRalston2`](@ref)                                      | 2      | 2     |
| [`TableauRalston3`](@ref)                                      | 3      | 3     |
| [`TableauRunge`](@ref), [`TableauRunge2`](@ref)                | 2      | 2     |
| [`TableauKutta`](@ref), [`TableauKutta3`](@ref)                | 3      | 3     |
| [`TableauRK416`](@ref), [`TableauRK4`](@ref)                   | 4      | 4     |
| [`TableauRK438`](@ref)                                         | 4      | 4     |
| [`TableauSSPRK3`](@ref)                                        | 3      | 3     |
| **Diagonally Implicit Methods**                                |        |       |
| [`TableauCrankNicolson`](@ref)                                 | 2      | 2     |
| [`TableauKraaijevangerSpijker`](@ref)                          | 2      | 2     |
| [`TableauQinZhang`](@ref)                                      | 2      | 2     |
| [`TableauCrouzeix`](@ref)                                      | 2      | 3     |
| **Fully Implicit Methods**                                     |        |       |
| [`TableauImplicitEuler`](@ref), [`TableauBackwardEuler`](@ref) | 1      | 1     |
| [`TableauImplicitMidpoint`](@ref)                              | 2      | 2     |
| [`TableauSRK3`](@ref)                                          | 3      | 4     |


## Gauß, Lobatto and Radau Methods

The coefficients of the Gauß, Lobatto and Radau methods are computed on-the-fly by the following constructors:

| Function                                    | Method                      | Order   |
|:--------------------------------------------|:----------------------------|:--------|
| [`TableauGLRK(s, T=Float64)`](@ref)         | Gauß-Legendre with s stages | $2s$    |
| [`TableauLobattoIIIA(s, T=Float64)`](@ref)  | Lobatto IIIA with s stages  | $2s-2$  |
| [`TableauLobattoIIIB(s, T=Float64)`](@ref)  | Lobatto IIIB with s stages  | $2s-2$  |
| [`TableauLobattoIIIC(s, T=Float64)`](@ref)  | Lobatto IIIC with s stages  | $2s-2$  |
| [`TableauLobattoIIIC̄(s, T=Float64)`](@ref)  | Lobatto IIIC̄ with s stages  | $2s-2$  |
| [`TableauLobattoIIID(s, T=Float64)`](@ref)  | Lobatto IIID with s stages  | $2s-2$  |
| [`TableauLobattoIIIE(s, T=Float64)`](@ref)  | Lobatto IIIE with s stages  | $2s-2$  |
| [`TableauLobattoIIIF(s, T=Float64)`](@ref)  | Lobatto IIIF with s stages  | $2s$    |
| [`TableauLobattoIIIG(s, T=Float64)`](@ref)  | Lobatto IIIG with s stages  | $2s$    |
| [`TableauRadauIA(s, T=Float64)`](@ref)      | Radau IA with s stages      | $2s-1$  |
| [`TableauRadauIIA(s, T=Float64)`](@ref)     | Radau IIA with s stages     | $2s-1$  |

The first argument `s` refers to the number of stages ($s \ge 1$ for Gauß and $s \ge 2$ for all other methods).
The second argument specifies the number type of the coefficients. Internally, all coefficients are computed using `BigFloat` and then converted to the requested number type, defaulting to `Float64`.
