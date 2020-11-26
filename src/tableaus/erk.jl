
"Tableau for explicit Euler method"
function TableauExplicitEuler(T=Float64)
    a = zeros(BigFloat, 1, 1)
    b = ones(BigFloat, 1)
    c = zeros(BigFloat, 1)
    o = 1

    Tableau{T}(:explicit_euler, o, a, b, c)
end

"Tableau for explicit midpoint method"
function TableauExplicitMidpoint(T=Float64)
    a = @big [[ 0.0  0.0 ]
              [ 0.5  0.0 ]]
    b = @big  [ 0.0, 1.0 ]
    c = @big  [ 0.0, 0.5 ]
    o = 2

    Tableau{T}(:explicit_midpoint, o, a, b, c)
end

"Tableau for Heun's method"
function TableauHeun(T=Float64)
    a = @big [[ 0.0  0.0 ]
              [ 1.0  0.0 ]]
    b = @big  [ 0.5, 0.5 ]
    c = @big  [ 0.0, 1.0 ]
    o = 2

    Tableau{T}(:heun, o, a, b, c)
end

"Tableau for Runge's method"
function TableauRunge(T=Float64)
    a = @big [[ 0.0  0.0 ]
              [ 1.0  0.0 ]]
    b = @big  [ 0.5, 0.5 ]
    c = @big  [ 0.0, 1.0 ]
    o = 2

    Tableau{T}(:runge, o, a, b, c)
end

"Tableau for Kutta's method of order three"
function TableauKutta(T=Float64)
    a = @big [[ 0.0  0.0  0.0 ]
              [ 0.5  0.0  0.0 ]
              [-1.0  2.0  0.0 ]]
    b = @big  [ 1/6, 4/6, 1/6 ]
    c = @big  [   0, 1/2,   1 ]
    o = 3

    Tableau{T}(:kutta, o, a, b, c)
end

"Tableau for explicit Runge-Kutta method of order four (1/6 rule)"
function TableauERK416(T=Float64)
    a = @big [[ 0.0  0.0  0.0  0.0 ]
              [ 0.5  0.0  0.0  0.0 ]
              [ 0.0  0.5  0.0  0.0 ]
              [ 0.0  0.0  1.0  0.0 ]]
    b = @big  [ 1/6, 1/3, 1/3, 1/6 ]
    c = @big  [ 0.0, 0.5, 0.5, 1.0 ]
    o = 4

    Tableau{T}(:erk4, o, a, b, c)
end

"Tableau for explicit Runge-Kutta method of order four (3/8 rule)"
function TableauERK438(T=Float64)
    a = @big [[ 0.0  0.0  0.0  0.0 ]
              [ 1/3  0.0  0.0  0.0 ]
              [-1/3  1.0  0.0  0.0 ]
              [ 1.0 -1.0  1.0  0.0 ]]
    b = @big  [ 1/8, 3/8, 3/8, 1/8 ]
    c = @big  [ 0.0, 1/3, 2/3, 1.0 ]
    o = 4

    Tableau{T}(:erk438, o, a, b, c)
end
