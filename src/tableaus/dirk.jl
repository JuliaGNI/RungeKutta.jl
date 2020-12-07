
"Tableau of Crank-Nicolson two-stage, 2nd order method"
function TableauCrankNicolson(::Type{T}) where {T}
    a = @big [[ 0     0   ]
              [ 1/2   1/2 ]]
    b = @big  [ 1/2,  1/2 ]
    c = @big  [ 0,    1   ]
    o = 2

    Tableau{T}(:cranknicolson, o, a, b, c)
end

TableauCrankNicolson() = TableauCrankNicolson(Float64)

"Tableau of Kraaijevanger and Spijker's two-stage, 2nd order method"
function TableauKraaijevangerSpijker(::Type{T}) where {T}
    a = @big [[ 1/2   0   ]
              [-1/2   2   ]]
    b = @big  [-1/2,  3/2 ]
    c = @big  [ 1/2,  3/2 ]
    o = 2

    Tableau{T}(:kraaijevangerspijker, o, a, b, c)
end

TableauKraaijevangerSpijker() = TableauKraaijevangerSpijker(Float64)

"Tableau of Qin and Zhang's symplectic two-stage, 2nd order method"
function TableauQinZhang(::Type{T}) where {T}
    a = @big [[ 1/4   0   ]
              [ 1/2   1/4 ]]
    b = @big  [ 1/2,  1/2 ]
    c = @big  [ 1/4,  3/4 ]
    o = 2

    Tableau{T}(:qinzhang, o, a, b, c)
end

TableauQinZhang() = TableauQinZhang(Float64)

"Tableau of Crouzeix's two-stage, 3rd order method"
function TableauCrouzeix(::Type{T}) where {T}
    fac = @big 1/2 / √3
    a = @big [[ 1/2+fac  0       ]
              [  -2*fac  1/2+fac ]]
    b = @big  [ 1/2,     1/2     ]
    c = @big  [ 1/2+fac, 1/2-fac ]
    o = 3

    Tableau{T}(:crouzeix, o, a, b, c)
end

TableauCrouzeix() = TableauCrouzeix(Float64)
