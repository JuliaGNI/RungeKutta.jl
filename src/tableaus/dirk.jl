
"Tableau of Crank-Nicolson two-stage, 2nd order method"
    a = @big [[ 0      0    ]
              [ 1//2   1//2 ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 0,     1    ]
function TableauCrankNicolson(::Type{T}) where {T}
    o = 2

    Tableau{T}(:cranknicolson, o, a, b, c)
end

TableauCrankNicolson() = TableauCrankNicolson(Float64)

"Tableau of Kraaijevanger and Spijker's two-stage, 2nd order method"
    a = @big [[ 1//2   0    ]
              [-1//2   2    ]]
    b = @big  [-1//2,  3//2 ]
    c = @big  [ 1//2,  3//2 ]
function TableauKraaijevangerSpijker(::Type{T}) where {T}
    o = 2

    Tableau{T}(:kraaijevangerspijker, o, a, b, c)
end

TableauKraaijevangerSpijker() = TableauKraaijevangerSpijker(Float64)

"Tableau of Qin and Zhang's symplectic two-stage, 2nd order method"
    a = @big [[ 1//4   0    ]
              [ 1//2   1//4 ]]
    b = @big  [ 1//2,  1//2 ]
    c = @big  [ 1//4,  3//4 ]
function TableauQinZhang(::Type{T}) where {T}
    o = 2

    Tableau{T}(:qinzhang, o, a, b, c)
end

TableauQinZhang() = TableauQinZhang(Float64)

"Tableau of Crouzeix's two-stage, 3rd order method"
    fac = @big 0.5/√3
    a = @big [[ 0.5+fac  0       ]
              [  -2*fac  0.5+fac ]]
    b = @big  [ 0.5,     0.5     ]
    c = @big  [ 0.5+fac, 0.5-fac ]
function TableauCrouzeix(::Type{T}) where {T}
    o = 3

    Tableau{T}(:crouzeix, o, a, b, c)
end

TableauCrouzeix() = TableauCrouzeix(Float64)
