
"Implicit Euler"
function TableauImplicitEuler()
    a = ones(Float64, 1, 1)
    b = [1.0]
    c = [1.0]
    o = 1

    Tableau(:implicit_euler, o, a, b, c)
end

"Implicit Midpoint"
function TableauImplicitMidpoint()
    a = 0.5*ones(Float64, 1, 1)
    b = [1.0]
    c = [0.5]
    o = 2

    Tableau(:implicit_midpoint, o, a, b, c)
end

"Symmetric Runge-Kutta tableau with three stages."
function TableauSRK3()
    a = [
         [5/36         2/9        5/36-√15/10]
         [5/36         2/9        5/36       ]
         [5/36+√15/10  2/9        5/36       ]
        ]
    b = [5/18,        4/9,       5/18       ]
    c = [1/2-√15/10,  1/2,       1/2+√15/10 ]
    o = 4

    Tableau{T}(:SRK3, o, a, b, c)
end
