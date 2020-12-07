
"Tableau of one-stage, 1st order implicit (backward) Euler method"
function TableauImplicitEuler(::Type{T}) where {T}
    a = ones(BigFloat, 1, 1)
    b = ones(BigFloat, 1)
    c = ones(BigFloat, 1)
    o = 1

    Tableau{T}(:implicit_euler, o, a, b, c)
end

TableauImplicitEuler() = TableauImplicitEuler(Float64)

"Alias for [`TableauImplicitEuler`](@ref)"
TableauBackwardEuler = TableauImplicitEuler

"Tableau of two-stage, 2nd order implicit midpoint method"
function TableauImplicitMidpoint(::Type{T}) where {T}
    a = ones(BigFloat, 1, 1) ./ 2
    b = ones(BigFloat, 1)
    c = ones(BigFloat, 1) ./ 2
    o = 2

    Tableau{T}(:implicit_midpoint, o, a, b, c)
end

TableauImplicitMidpoint() = TableauImplicitMidpoint(Float64)

"Tableau of symmetric and symplectic three-stage, 4th order Runge-Kutta method"
function TableauSRK3(::Type{T}) where {T}
    a = @big [[ 5/36         2/9        5/36-√15/10 ]
              [ 5/36         2/9        5/36        ]
              [ 5/36+√15/10  2/9        5/36        ]]
    b = @big  [ 5/18,        4/9,       5/18        ]
    c = @big  [ 1/2-√15/10,  1/2,       1/2+√15/10  ]
    o = 4

    Tableau{T}(:SRK3, o, a, b, c)
end

TableauSRK3() = TableauSRK3(Float64)
