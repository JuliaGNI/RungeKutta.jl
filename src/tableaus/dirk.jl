
function TableauCrouzeix(T=Float64)
    fac = @big 0.5/√3
    a = @big [[ 0.5+fac  0       ]
              [  -2*fac  0.5+fac ]]
    b = @big  [ 0.5,     0.5     ]
    c = @big  [ 0.5+fac, 0.5-fac ]
    o = 3

    Tableau{T}(:crouzeix, o, a, b, c)
end
