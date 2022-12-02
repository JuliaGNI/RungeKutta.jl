```@meta
CurrentModule = RungeKutta
```

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
