
GeometricBase.reference(::Val{:CrankNicolson}) = """
Reference:

    J. Crank and P. Nicolson.
    A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type.
    Mathematical Proceedings of the Cambridge Philosophical Society, Volume 43, Issue 1, Pages 50-67, 1947.
    doi: 10.1017/S0305004100023197
"""

"""
Tableau of Crank-Nicolson two-stage, 2nd order method

```julia
TableauCrankNicolson(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:CrankNicolson)))
"""
function TableauCrankNicolson(::Type{T}=Float64) where {T}
    a = @big [[ 0     0   ]
              [ 1/2   1/2 ]]
    b = @big  [ 1/2,  1/2 ]
    c = @big  [ 0,    1   ]
    o = 2

    Tableau{T}(:CrankNicolson, o, a, b, c)
end


GeometricBase.reference(::Val{:KraaijevangerSpijker}) = """
Reference:

    J. F. B. M. Kraaijevanger and M. N. Spijker.
    Algebraic stability and error propagation in Runge-Kutta methods.
    Applied Numerical Mathematics, Volume 5, Issues 1-2, Pages 71-87, 1989.
    doi: 10.1016/0168-9274(89)90025-1
"""

"""
Tableau of Kraaijevanger and Spijker's two-stage, 2nd order method

```julia
TableauKraaijevangerSpijker(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:KraaijevangerSpijker)))
"""
function TableauKraaijevangerSpijker(::Type{T}=Float64) where {T}
    a = @big [[ 1/2   0   ]
              [-1/2   2   ]]
    b = @big  [-1/2,  3/2 ]
    c = @big  [ 1/2,  3/2 ]
    o = 2

    Tableau{T}(:KraaijevangerSpijker, o, a, b, c)
end


GeometricBase.reference(::Val{:QinZhang}) = """
Reference:

    M.-Z. Qin and M.-Q. Zhang.
    Symplectic Runge-Kutta algorithms for Hamilton systems.
    Journal of Computational Mathematics, Supplementary Issue, Pages 205-215, 1992.
"""

"""
Tableau of Qin and Zhang's symplectic two-stage, 2nd order method

```julia
TableauQinZhang(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:QinZhang)))
"""
function TableauQinZhang(::Type{T}=Float64) where {T}
    a = @big [[ 1/4   0   ]
              [ 1/2   1/4 ]]
    b = @big  [ 1/2,  1/2 ]
    c = @big  [ 1/4,  3/4 ]
    o = 2

    Tableau{T}(:QinZhang, o, a, b, c)
end


GeometricBase.reference(::Val{:Crouzeix}) = """
Reference:

    M.Crouzeix.
    Sur L'approximation des équations différentielles opérationelles linéaires par des méthodes de Runge-Kutta.
    Thesis. Université de Paris, 1975.
"""

"""
Tableau of Crouzeix's two-stage, 3rd order method

```julia
TableauCrouzeix(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:Crouzeix)))
"""
function TableauCrouzeix(::Type{T}=Float64) where {T}
    fac = @big 1/2 / √3
    a = @big [[ 1/2+fac  0       ]
              [  -2*fac  1/2+fac ]]
    b = @big  [ 1/2,     1/2     ]
    c = @big  [ 1/2+fac, 1/2-fac ]
    o = 3

    Tableau{T}(:Crouzeix, o, a, b, c)
end
