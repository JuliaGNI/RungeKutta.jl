
"""
Tableau of one-stage, 1st order explicit (forward) Euler method

Reference:

    Leonhard Euler
    Institutiones calculi differentialis cum eius vsu in analysi finitorum ac doctrina serierum.
    Imp. Acad. Imper. Scient. Petropolitanae, Opera Omnia, Vol.X, [I.6], 1755.
    In: Opera Omnia, 1st Series, Volume 11, Institutiones Calculi Integralis. Teubner, Leipzig, Pages 424-434, 1913.
    Sectio secunda. Caput VII. De integratione aequationum differentialium per approximationem. Problema 85.

"""
function TableauExplicitEuler(::Type{T}=Float64) where {T}
    a = zeros(BigFloat, 1, 1)
    b = ones(BigFloat, 1)
    c = zeros(BigFloat, 1)
    o = 1

    Tableau{T}(:explicit_euler, o, a, b, c)
end

"Alias for [`TableauExplicitEuler`](@ref)"
TableauForwardEuler = TableauExplicitEuler

"""
Tableau of explicit two-stage, 2nd order midpoint method

Reference:

    Carl Runge
    Über die numerische Auflösung von Differentialgleichungen.
    Mathematische Annalen, Volume 46, Pages 167-178, 1895.
    doi: 10.1007/BF01446807.
    Equation (2)

"""
function TableauExplicitMidpoint(::Type{T}=Float64) where {T}
    a = @big [[ 0     0    ]
              [ 1//2  0    ]]
    b = @big  [ 0,    1    ]
    c = @big  [ 0,    1//2 ]
    o = 2

    Tableau{T}(:explicit_midpoint, o, a, b, c)
end

"""
Tableau of Heun's two-stage, 2nd order method

Reference:

    Karl Heun.
    Neue Methoden zur approximativen Integration der Differentialgleichungen einer unabhängigen Veränderlichen.
    Zeitschrift für Mathematik und Physik, Volume 45, Pages 23-38, 1900.
    Algorithm II.

"""
function TableauHeun2(::Type{T}=Float64) where {T}
    a = @big [[ 0      0    ]
              [ 1      0    ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 0,     1    ]
    o = 2

    Tableau{T}(:heun2, o, a, b, c)
end

"""
Tableau of Heun's three-stage, 3rd order method

Reference:

    Karl Heun.
    Neue Methoden zur approximativen Integration der Differentialgleichungen einer unabhängigen Veränderlichen.
    Zeitschrift für Mathematik und Physik, Volume 45, Pages 23-38, 1900.
    Algorithm VI.

"""
function TableauHeun3(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//3   0      0    ]
              [ 0      2//3   0    ]]
    b = @big  [ 1//4,  0,     3//4 ]
    c = @big  [ 0,     1//3,  2//3 ]
    o = 3

    Tableau{T}(:heun3, o, a, b, c)
end

"""
Tableau of Ralston's two-stage, 2nd order method

Reference:

    Anthony Ralston.
    Runge-Kutta Methods with Minimum Error Bounds.
    Mathematics of Computation, Volume 16, Pages 431-437, 1962.
    doi: 10.1090/S0025-5718-1962-0150954-0.
    Equation (3.5)

"""
function TableauRalston2(::Type{T}=Float64) where {T}
    a = @big [[ 0      0    ]
              [ 2//3   0    ]]
    b = @big  [ 1//4,  3//4 ]
    c = @big  [ 0,     2//3 ]
    o = 2

    Tableau{T}(:ralston2, o, a, b, c)
end

"""
Tableau of Ralston's three-stage, 3rd order method

Reference:

    Anthony Ralston.
    Runge-Kutta Methods with Minimum Error Bounds.
    Mathematics of Computation, Volume 16, Pages 431-437, 1962.
    doi: 10.1090/S0025-5718-1962-0150954-0.
    Equation (4.10)

"""
function TableauRalston3(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//2   0      0    ]
              [ 0      3//4   0    ]]
    b = @big  [ 2//9,  3//9,  4//9 ]
    c = @big  [ 0,     1//2,  3//4 ]
    o = 3

    Tableau{T}(:ralston3, o, a, b, c)
end

"""
Tableau of Runge's two-stage, 2nd order method

Reference:

    Carl Runge
    Über die numerische Auflösung von Differentialgleichungen.
    Mathematische Annalen, Volume 46, Pages 167-178, 1895.
    doi: 10.1007/BF01446807.
    Equation (3)

"""
function TableauRunge(::Type{T}=Float64) where {T}
    a = @big [[ 0      0    ]
              [ 1      0    ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 0,     1    ]
    o = 2

    Tableau{T}(:runge, o, a, b, c)
end

"Alias for [`TableauRunge`](@ref)"
TableauRunge2 = TableauRunge

"""
Tableau of Kutta's three-stage, 3rd order method

Reference:

    Wilhelm Kutta
    Beitrag zur Näherungsweisen Integration totaler Differentialgleichungen
    Zeitschrift für Mathematik und Physik, Volume 46, Pages 435–453, 1901.
    Page 440

"""
function TableauKutta(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//2   0      0    ]
              [-1      2      0    ]]
    b = @big  [ 1//6,  4//6,  1//6 ]
    c = @big  [ 0,     1//2,  1    ]
    o = 3

    Tableau{T}(:kutta, o, a, b, c)
end

"Alias for [`TableauKutta`](@ref)"
TableauKutta3 = TableauKutta

"""
Tableau of explicit Runge-Kutta method of order four (1/6 rule)

Reference:

    Wilhelm Kutta
    Beitrag zur Näherungsweisen Integration totaler Differentialgleichungen
    Zeitschrift für Mathematik und Physik, Volume 46, Pages 435–453, 1901.
    Page 443

"""
function TableauRK416(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0      0    ]
              [ 1//2   0      0      0    ]
              [ 0      1//2   0      0    ]
              [ 0      0      1      0    ]]
    b = @big  [ 1//6,  1//3,  1//3,  1//6 ]
    c = @big  [ 0,     1//2,  1//2,  1    ]
    o = 4

    Tableau{T}(:erk4, o, a, b, c)
end

"Alias for [`TableauRK416`](@ref)"
TableauRK4 = TableauRK416

"""
Tableau of explicit Runge-Kutta method of order four (3/8 rule)

Reference:

    Wilhelm Kutta
    Beitrag zur Näherungsweisen Integration totaler Differentialgleichungen
    Zeitschrift für Mathematik und Physik, Volume 46, Pages 435–453, 1901.
    Page 441

"""
function TableauRK438(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0      0    ]
              [ 1//3   0      0      0    ]
              [-1//3   1      0      0    ]
              [ 1     -1      1      0    ]]
    b = @big  [ 1//8,  3//8,  3//8,  1//8 ]
    c = @big  [ 0,     1//3,  2//3,  1    ]
    o = 4

    Tableau{T}(:erk438, o, a, b, c)
end

"""
Tableau of 2rd order Strong Stability Preserving method with two stages and CFL ≤ 1

Alias for [`TableauHeun2`](@ref)

Reference:

    Chi-Wang Shu, Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics, Volume 77, Issue 2, Pages 439-471, 1988.
    doi: 10.1016/0021-9991(88)90177-5.
    Equation (2.16)

"""
TableauSSPRK2 = TableauHeun2

"""
Tableau of 3rd order Strong Stability Preserving method with three stages and CFL ≤ 1

Reference:

    Chi-Wang Shu, Stanley Osher.
    Efficient implementation of essentially non-oscillatory shock-capturing schemes.
    Journal of Computational Physics, Volume 77, Issue 2, Pages 439-471, 1988.
    doi: 10.1016/0021-9991(88)90177-5.
    Equation (2.18)

"""
function TableauSSPRK3(::Type{T}=Float64) where {T}
    a = @big [[ 0      0      0    ]
              [ 1      0      0    ]
              [ 1//4   1//4   0    ]]
    b = @big  [ 1//6,  1//6,  4//6 ]
    c = @big  [ 0,     1,     1//2 ]
    o = 3

    Tableau{T}(:ssprk3, o, a, b, c)
end
