
reference(::Val{:ImplicitEuler}) = """
Reference:

    Augustin-Louis Cauchy.
    Équations différentielles ordinaires. Cours inédit (fragment). Douzième leçon.
    Ed. Christian Gilain, Etudes Vivantes, 1981.
    Page 102, Equation (5), Θ=1.
"""

"""
Tableau of one-stage, 1st order implicit (backward) Euler method

```julia
TableauImplicitEuler(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:ImplicitEuler)))
"""
function TableauImplicitEuler(::Type{T}=Float64) where {T}
    a = ones(BigFloat, 1, 1)
    b = ones(BigFloat, 1)
    c = ones(BigFloat, 1)
    o = 1

    Tableau{T}(:ImplicitEuler, o, a, b, c; R∞=0)
end

"Alias for [`TableauImplicitEuler`](@ref)"
const TableauBackwardEuler = TableauImplicitEuler
reference(::Val{:BackwardEuler}) = reference(Val(:ImplicitEuler))


reference(::Val{:ImplicitMidpoint}) = """
Reference:

    Augustin-Louis Cauchy.
    Équations différentielles ordinaires. Cours inédit (fragment). Douzième leçon.
    Ed. Christian Gilain, Etudes Vivantes, 1981.
    Page 102, Equation (5), Θ=1/2.
"""

"""
Tableau of two-stage, 2nd order implicit midpoint method

```julia
TableauImplicitMidpoint(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:ImplicitMidpoint)))
"""
function TableauImplicitMidpoint(::Type{T}=Float64) where {T}
    a = ones(BigFloat, 1, 1) ./ 2
    b = ones(BigFloat, 1)
    c = ones(BigFloat, 1) ./ 2
    o = 2

    Tableau{T}(:ImplicitMidpoint, o, a, b, c; R∞=-1)
end


reference(::Val{:IRK3}) = """
Reference:

    Ernst Hairer and Gerhard Wanner.
    Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems.
    Springer, 1996.
    Section IV.5, W-transformation.
    The s=2, γ=1/2 member of the family of algebraically stable methods of order 2s-1,
    obtained from the Gauss method by X = Xₛ + γ eₛ eₛᵀ.
"""

"""
Tableau of two-stage, 3rd order fully implicit Runge-Kutta method

```julia
TableauIRK3(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

The method uses the two-point Gauss-Legendre nodes and weights, but its coefficient matrix is
that of [`TableauGauss`](@ref) plus the rank-one perturbation `¼ (e₁-e₂) (e₁-e₂)ᵀ`. It is
therefore of order three rather than four and neither symmetric nor symplectic, but it is
A-stable and algebraically stable with `R(∞) = -1/2`.

$(reference(Val(:IRK3)))
"""
function TableauIRK3(::Type{T}=Float64) where {T}
    a = @big [[ 1/2        -√3/6     ]
              [+√3/6        1/2      ]]
    b = @big  [ 1/2,        1/2      ]
    c = @big  [ 1/2-√3/6,   1/2+√3/6 ]
    o = 3

    Tableau{T}(:IRK3, o, a, b, c; R∞=-1//2)
end


reference(::Val{:SRK3}) = """
Reference:

    Shan Zhao and Guo-Wei Wei.
    A unified discontinuous Galerkin framework for time integration.
    Mathematical Methods in the Applied Sciences, Volume 37, Issue 7, Pages 1042-1071, 2014.
    doi: 10.1002/mma.2863.
"""

"""
Tableau of symmetric and symplectic three-stage, 4th order Runge-Kutta method

```julia
TableauSRK3(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:SRK3)))
"""
function TableauSRK3(::Type{T}=Float64) where {T}
    a = @big [[ 5/36         2/9        5/36-√15/10 ]
              [ 5/36         2/9        5/36        ]
              [ 5/36+√15/10  2/9        5/36        ]]
    b = @big  [ 5/18,        4/9,       5/18        ]
    c = @big  [ 1/2-√15/10,  1/2,       1/2+√15/10  ]
    o = 4

    Tableau{T}(:SRK3, o, a, b, c; R∞=-1)
end
