
RungeKutta.reference(::Val{:ImplicitEuler}) = """
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
TableauBackwardEuler = TableauImplicitEuler


RungeKutta.reference(::Val{:ImplicitMidpoint}) = """
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


RungeKutta.reference(::Val{:SRK3}) = """
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
