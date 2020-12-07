
"Tableau of one-stage, 1st order explicit (forward) Euler method"
function TableauExplicitEuler(::Type{T}) where {T}
    a = zeros(BigFloat, 1, 1)
    b = ones(BigFloat, 1)
    c = zeros(BigFloat, 1)
    o = 1

    Tableau{T}(:explicit_euler, o, a, b, c)
end

TableauExplicitEuler() = TableauExplicitEuler(Float64)

"Alias for [`TableauExplicitEuler`](@ref)"
TableauForwardEuler = TableauExplicitEuler

"Tableau of explicit two-stage, 2nd order midpoint method"
function TableauExplicitMidpoint(::Type{T}) where {T}
    a = @big [[ 0     0    ]
              [ 1//2  0    ]]
    b = @big  [ 0,    1    ]
    c = @big  [ 0,    1//2 ]
    o = 2

    Tableau{T}(:explicit_midpoint, o, a, b, c)
end

TableauExplicitMidpoint() = TableauExplicitMidpoint(Float64)

"Tableau of Heun's two-stage, 2nd order method"
function TableauHeun2(::Type{T}) where {T}
    a = @big [[ 0      0    ]
              [ 1      0    ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 0,     1    ]
    o = 2

    Tableau{T}(:heun2, o, a, b, c)
end

TableauHeun2() = TableauHeun2(Float64)

"Tableau of Heun's three-stage, 3rd order method"
function TableauHeun3(::Type{T}) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//3   0      0    ]
              [ 0      2//3   0    ]]
    b = @big  [ 1//4,  0,     3//4 ]
    c = @big  [ 0,     1//3,  2//3 ]
    o = 3

    Tableau{T}(:heun3, o, a, b, c)
end

TableauHeun3() = TableauHeun3(Float64)

"Tableau of Ralston's two-stage, 2nd order method"
function TableauRalston2(::Type{T}) where {T}
    a = @big [[ 0      0    ]
              [ 2//3   0    ]]
    b = @big  [ 1//4,  3//4 ]
    c = @big  [ 0,     2//3 ]
    o = 2

    Tableau{T}(:ralston2, o, a, b, c)
end

TableauRalston2() = TableauRalston2(Float64)

"Tableau of Ralston's three-stage, 3rd order method"
function TableauRalston3(::Type{T}) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//2   0      0    ]
              [ 0      3//4   0    ]]
    b = @big  [ 2//9,  3//9,  4//9 ]
    c = @big  [ 0,     1//2,  3//4 ]
    o = 3

    Tableau{T}(:ralston3, o, a, b, c)
end

TableauRalston3() = TableauRalston3(Float64)

"Tableau of Runge's two-stage, 2nd order method"
function TableauRunge(::Type{T}) where {T}
    a = @big [[ 0      0    ]
              [ 1      0    ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 0,     1    ]
    o = 2

    Tableau{T}(:runge, o, a, b, c)
end

TableauRunge() = TableauRunge(Float64)

"Alias for [`TableauRunge`](@ref)"
TableauRunge2 = TableauRunge

"Tableau of Kutta's three-stage, 3rd order method"
function TableauKutta(::Type{T}) where {T}
    a = @big [[ 0      0      0    ]
              [ 1//2   0      0    ]
              [-1      2      0    ]]
    b = @big  [ 1//6,  4//6,  1//6 ]
    c = @big  [ 0,     1//2,  1    ]
    o = 3

    Tableau{T}(:kutta, o, a, b, c)
end

TableauKutta() = TableauKutta(Float64)

"Alias for [`TableauKutta`](@ref)"
TableauKutta3 = TableauKutta

"Tableau of explicit Runge-Kutta method of order four (1/6 rule)"
function TableauRK416(::Type{T}) where {T}
    a = @big [[ 0      0      0      0    ]
              [ 1//2   0      0      0    ]
              [ 0      1//2   0      0    ]
              [ 0      0      1      0    ]]
    b = @big  [ 1//6,  1//3,  1//3,  1//6 ]
    c = @big  [ 0,     1//2,  1//2,  1    ]
    o = 4

    Tableau{T}(:erk4, o, a, b, c)
end

TableauRK416() = TableauRK416(Float64)

"Alias for [`TableauRK416`](@ref)"
TableauRK4 = TableauRK416

"Tableau of explicit Runge-Kutta method of order four (3/8 rule)"
function TableauRK438(::Type{T}) where {T}
    a = @big [[ 0      0      0      0    ]
              [ 1//3   0      0      0    ]
              [-1//3   1      0      0    ]
              [ 1     -1      1      0    ]]
    b = @big  [ 1//8,  3//8,  3//8,  1//8 ]
    c = @big  [ 0,     1//3,  2//3,  1    ]
    o = 4

    Tableau{T}(:erk438, o, a, b, c)
end

TableauRK438() = TableauRK438(Float64)

"Tableau of 3rd order Strong Stability Preserving method with three stages"
function TableauSSPRK3(::Type{T}) where {T}
    a = @big [[ 0      0      0    ]
              [ 1      0      0    ]
              [ 1//4   1//4   0    ]]
    b = @big  [ 1//6,  1//6,  4//6 ]
    c = @big  [ 0,     1,     1//2 ]
    o = 3

    Tableau{T}(:ssprk3, o, a, b, c)
end

TableauSSPRK3() = TableauSSPRK3(Float64)
