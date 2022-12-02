
RungeKutta.reference(::Val{:ExplicitEuler}) = """
Reference:

    Leonhard Euler.
    Institutiones calculi differentialis cum eius vsu in analysi finitorum ac doctrina serierum.
    Imp. Acad. Imper. Scient. Petropolitanae, Opera Omnia, Vol.X, [I.6], 1755.
    In: Opera Omnia, 1st Series, Volume 11, Institutiones Calculi Integralis. Teubner, Leipzig, Pages 424-434, 1913.
    Sectio secunda. Caput VII. De integratione aequationum differentialium per approximationem. Problema 85.
"""

"""
Tableau of one-stage, 1st order explicit (forward) Euler method

```julia
TableauExplicitEuler(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:ExplicitEuler)))
"""
function TableauExplicitEuler(::Type{T}=Float64) where {T}
    a = zeros(BigFloat, 1, 1)
    b = ones(BigFloat, 1)
    c = zeros(BigFloat, 1)
    o = 1

    Tableau{T}(:ExplicitEuler, o, a, b, c; R∞=Inf)
end

"Alias for [`TableauExplicitEuler`](@ref)"
TableauForwardEuler = TableauExplicitEuler


RungeKutta.reference(::Val{:ExplicitMidpoint}) = """
Reference:

    Carl Runge.
    Über die numerische Auflösung von Differentialgleichungen.
    Mathematische Annalen, Volume 46, Pages 167-178, 1895.
    doi: 10.1007/BF01446807.
    Equation (2)
"""

"""
Tableau of explicit two-stage, 2nd order midpoint method

```julia
TableauExplicitMidpoint(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:ExplicitMidpoint)))
"""
function TableauExplicitMidpoint(::Type{T}=Float64) where {T}
    a = @big [[ 0     0    ]
              [ 1//2  0    ]]
    b = @big  [ 0,    1    ]
    c = @big  [ 0,    1//2 ]
    o = 2

    Tableau{T}(:ExplicitMidpoint, o, a, b, c)
end


RungeKutta.reference(::Val{:Heun2}) = """
Reference:

    Karl Heun.
    Neue Methoden zur approximativen Integration der Differentialgleichungen einer unabhängigen Veränderlichen.
    Zeitschrift für Mathematik und Physik, Volume 45, Pages 23-38, 1900.
    Algorithm II.
"""

"""
Tableau of Heun's two-stage, 2nd order method

```julia
TableauHeun2(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:Heun2)))
"""
function TableauHeun2(::Type{T}=Float64) where {T}
    a = @big [[ 0      0    ]
              [ 1      0    ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 0,     1    ]
    o = 2

    Tableau{T}(:Heun2, o, a, b, c)
end


RungeKutta.reference(::Val{:Heun3}) = """
Reference:

    Karl Heun.
    Neue Methoden zur approximativen Integration der Differentialgleichungen einer unabhängigen Veränderlichen.
    Zeitschrift für Mathematik und Physik, Volume 45, Pages 23-38, 1900.
    Algorithm VI.
"""

"""
Tableau of Heun's three-stage, 3rd order method

```julia
TableauHeun3(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:Heun3)))
"""
function TableauHeun3(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//3   0      0    ]
              [ 0      2//3   0    ]]
    b = @big  [ 1//4,  0,     3//4 ]
    c = @big  [ 0,     1//3,  2//3 ]
    o = 3

    Tableau{T}(:Heun3, o, a, b, c)
end


RungeKutta.reference(::Val{:Ralston2}) = """
Reference:

    Anthony Ralston.
    Runge-Kutta Methods with Minimum Error Bounds.
    Mathematics of Computation, Volume 16, Pages 431-437, 1962.
    doi: 10.1090/S0025-5718-1962-0150954-0.
    Equation (3.5)
"""

"""
Tableau of Ralston's two-stage, 2nd order method

```julia
TableauRalston2(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:Ralston2)))
"""
function TableauRalston2(::Type{T}=Float64) where {T}
    a = @big [[ 0      0    ]
              [ 2//3   0    ]]
    b = @big  [ 1//4,  3//4 ]
    c = @big  [ 0,     2//3 ]
    o = 2

    Tableau{T}(:Ralston2, o, a, b, c)
end


RungeKutta.reference(::Val{:Ralston3}) = """
Reference:

    Anthony Ralston.
    Runge-Kutta Methods with Minimum Error Bounds.
    Mathematics of Computation, Volume 16, Pages 431-437, 1962.
    doi: 10.1090/S0025-5718-1962-0150954-0.
    Equation (4.10)
"""

"""
Tableau of Ralston's three-stage, 3rd order method

```julia
TableauRalston3(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:Ralston3)))
"""
function TableauRalston3(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//2   0      0    ]
              [ 0      3//4   0    ]]
    b = @big  [ 2//9,  3//9,  4//9 ]
    c = @big  [ 0,     1//2,  3//4 ]
    o = 3

    Tableau{T}(:Ralston3, o, a, b, c)
end


RungeKutta.reference(::Val{:Runge}) = """
Reference:

    Carl Runge
    Über die numerische Auflösung von Differentialgleichungen.
    Mathematische Annalen, Volume 46, Pages 167-178, 1895.
    doi: 10.1007/BF01446807.
    Equation (3)
"""

"""
Tableau of Runge's two-stage, 2nd order method

```julia
TableauRunge(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:Runge)))
"""
function TableauRunge(::Type{T}=Float64) where {T}
    a = @big [[ 0      0    ]
              [ 1      0    ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 0,     1    ]
    o = 2

    Tableau{T}(:Runge, o, a, b, c)
end

"Alias for [`TableauRunge`](@ref)"
TableauRunge2 = TableauRunge


RungeKutta.reference(::Val{:Kutta}) = """
Reference:

    Wilhelm Kutta
    Beitrag zur Näherungsweisen Integration totaler Differentialgleichungen
    Zeitschrift für Mathematik und Physik, Volume 46, Pages 435–453, 1901.
    Page 440
"""

"""
Tableau of Kutta's three-stage, 3rd order method

```julia
TableauKutta(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:Kutta)))
"""
function TableauKutta(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//2   0      0    ]
              [-1      2      0    ]]
    b = @big  [ 1//6,  4//6,  1//6 ]
    c = @big  [ 0,     1//2,  1    ]
    o = 3

    Tableau{T}(:Kutta, o, a, b, c)
end

"Alias for [`TableauKutta`](@ref)"
TableauKutta3 = TableauKutta


RungeKutta.reference(::Val{:RK416}) = """
Reference:

    Wilhelm Kutta
    Beitrag zur Näherungsweisen Integration totaler Differentialgleichungen
    Zeitschrift für Mathematik und Physik, Volume 46, Pages 435–453, 1901.
    Page 443
"""

"""
Tableau of explicit Runge-Kutta method of order four (1/6 rule)

```julia
TableauRK416(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:RK416)))
"""
function TableauRK416(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0      0    ]
              [ 1//2   0      0      0    ]
              [ 0      1//2   0      0    ]
              [ 0      0      1      0    ]]
    b = @big  [ 1//6,  1//3,  1//3,  1//6 ]
    c = @big  [ 0,     1//2,  1//2,  1    ]
    o = 4

    Tableau{T}(:RK416, o, a, b, c)
end

"Alias for [`TableauRK416`](@ref)"
TableauRK4 = TableauRK416


RungeKutta.reference(::Val{:RK438}) = """
Reference:

    Wilhelm Kutta
    Beitrag zur Näherungsweisen Integration totaler Differentialgleichungen
    Zeitschrift für Mathematik und Physik, Volume 46, Pages 435–453, 1901.
    Page 441
"""

"""
Tableau of explicit Runge-Kutta method of order four (3/8 rule)

```julia
TableauRK438(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:RK438)))
"""
function TableauRK438(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0      0    ]
              [ 1//3   0      0      0    ]
              [-1//3   1      0      0    ]
              [ 1     -1      1      0    ]]
    b = @big  [ 1//8,  3//8,  3//8,  1//8 ]
    c = @big  [ 0,     1//3,  2//3,  1    ]
    o = 4

    Tableau{T}(:RK438, o, a, b, c)
end


RungeKutta.reference(::Val{:SSPRK2}) = """
Reference:

    Chi-Wang Shu, Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics, Volume 77, Issue 2, Pages 439-471, 1988.
    doi: 10.1016/0021-9991(88)90177-5.
    Equation (2.16)
"""

"""
Tableau of 2rd order Strong Stability Preserving method with two stages and CFL ≤ 1

```julia
TableauSSPRK2(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

This is the same tableau as [`TableauHeun2`](@ref).

$(reference(Val(:SSPRK2)))
"""
function TableauSSPRK2(args...)
    tab = TableauHeun2(args...)
    Tableau(:SSPRK2, tab.o, tab.a, tab.b, tab.c)
end


RungeKutta.reference(::Val{:SSPRK3}) = """
Reference:

    Chi-Wang Shu, Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics, Volume 77, Issue 2, Pages 439-471, 1988.
    doi: 10.1016/0021-9991(88)90177-5.
    Equation (2.18)
"""

"""
Tableau of 3rd order Strong Stability Preserving method with three stages and CFL ≤ 1

```julia
TableauSSPRK3(::Type{T}=Float64) where {T}
```
The constructor takes one optional argument, that is the element type of the tableau.

$(reference(Val(:SSPRK3)))
"""
function TableauSSPRK3(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1      0      0    ]
              [ 1//4   1//4   0    ]]
    b = @big  [ 1//6,  1//6,  4//6 ]
    c = @big  [ 0,     1,     1//2 ]
    o = 3

    Tableau{T}(:SSPRK3, o, a, b, c)
end
