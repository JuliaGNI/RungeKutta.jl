module Tableaus

    import Polynomials
    import Polynomials: Polynomial

    using ..RungeKutta
    using ..RungeKutta: big, @big, _legendre, _shifted_legendre
    

    include("tableaus/dirk.jl")

    export TableauCrankNicolson,
           TableauCrouzeix,
           TableauKraaijevangerSpijker,
           TableauQinZhang

    include("tableaus/erk.jl")

    export TableauExplicitEuler, TableauForwardEuler,
           TableauExplicitMidpoint, 
           TableauHeun2, TableauRK21,
           TableauHeun3,
           TableauRalston2,
           TableauRalston3,
           TableauRunge, TableauRunge2, TableauRK22,
           TableauKutta, TableauKutta3, TableauRK32,
           TableauRK31,
           TableauRK416, TableauRK41, TableauRK4,
           TableauRK42,
           TableauRK438,
           TableauRK5,
           TableauSSPRK2,
           TableauSSPRK3

    include("tableaus/firk.jl")

    export TableauImplicitEuler, TableauBackwardEuler,
           TableauImplicitMidpoint,
           TableauSRK3

    include("tableaus/gauss.jl")
    
    export TableauGauss

    include("tableaus/lobatto.jl")

    export TableauLobattoIII,
           TableauLobattoIIIA,
           TableauLobattoIIIĀ,
           TableauLobattoIIIB,
           TableauLobattoIIIB̄,
           TableauLobattoIIIC,
           TableauLobattoIIIC̄,
           TableauLobattoIIID,
           TableauLobattoIIID̄,
           TableauLobattoIIIE,
           TableauLobattoIIIĒ,
           TableauLobattoIIIF,
           TableauLobattoIIIF̄,
           TableauLobattoIIIG,
           TableauLobattoIIIḠ
    export get_lobatto_nullvector

    include("tableaus/radau.jl")

    export TableauRadauIA, TableauRadauIIA,
           TableauRadauIB, TableauRadauIIB


    TableauList = (
       # explicit
       forward_euler         = TableauForwardEuler,
       explicit_euler        = TableauExplicitEuler,
       explicit_midpoint     = TableauExplicitMidpoint, 
       heun2                 = TableauHeun2,
       heun3                 = TableauHeun3,
       ralston2              = TableauRalston2,
       ralston3              = TableauRalston3,
       runge                 = TableauRunge,
       runge2                = TableauRunge2,
       kutta                 = TableauKutta,
       kutta3                = TableauKutta3,
       rk21                  = TableauRK21,
       rk22                  = TableauRK22,
       rk31                  = TableauRK31,
       rk32                  = TableauRK32,
       rk4                   = TableauRK4,
       rk41                  = TableauRK41,
       rk416                 = TableauRK416,
       rk42                  = TableauRK42,
       rk438                 = TableauRK438,
       rk5                   = TableauRK5,
       ssprk2                = TableauSSPRK2,
       ssprk3                = TableauSSPRK3,
       # dirk
       crank_nicolson        = TableauCrankNicolson,
       crouzeix              = TableauCrouzeix,
       kraaijevanger_spijker = TableauKraaijevangerSpijker,
       qin_zhang             = TableauQinZhang,
       # firk
       backward_euler        = TableauBackwardEuler,
       implicit_euler        = TableauImplicitEuler,
       implicit_midpoint     = TableauImplicitMidpoint,
       srk3                  = TableauSRK3,
    )

    export TableauList

end
