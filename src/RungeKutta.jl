module RungeKutta

    import GenericLinearAlgebra
    import Polynomials
    import Polynomials: Polynomial

    include("utils.jl")

    include("tableau.jl")

    export Tableau,
           isexplicit,
           isimplicit,
           isdiagnonallyimplicit,
           isfullyimplicit

    include("order_conditions.jl")

    export check_order_conditions_b,
           check_order_conditions_c,
           check_order_conditions_d,
           satisfies_simplifying_assumption_b,
           satisfies_simplifying_assumption_c,
           satisfies_simplifying_assumption_d,
           solve_simplifying_assumption_b,
           solve_simplifying_assumption_c,
           solve_simplifying_assumption_d

    include("symmetry.jl")
    
    export check_symmetry,
           issymmetric

    include("symplecticity.jl")
    
    export check_symplecticity,
           compute_symplecticity_error,
           issymplectic,
           get_symplectic_conjugate_coefficients,
           get_symplectic_conjugate_coefficients!,
           symplecticize

    include("tableaus/dirk.jl")

    export TableauKraaijevangerSpijker,
           TableauQinZhang,
           TableauCrouzeix

    include("tableaus/erk.jl")

    export TableauExplicitEuler, TableauForwardEuler,
           TableauExplicitMidpoint, 
           TableauHeun2,
           TableauHeun3,
           TableauRalston2,
           TableauRalston3,
           TableauRunge, TableauRunge2,
           TableauKutta, TableauKutta3,
           TableauRK416, TableauRK4,
           TableauRK438,
           TableauSSPRK2,
           TableauSSPRK3

    include("tableaus/firk.jl")

    export TableauImplicitEuler, TableauBackwardEuler,
           TableauImplicitMidpoint,
           TableauCrankNicolson,
           TableauSRK3

    include("tableaus/gauss.jl")
    
    export TableauGauss

    include("tableaus/lobatto.jl")

    export TableauLobattoIIIA, TableauLobattoIIIB, TableauLobattoIIIC, TableauLobattoIIIC̄,
           TableauLobattoIIID, TableauLobattoIIIE, TableauLobattoIIIF, TableauLobattoIIIG

    include("tableaus/radau.jl")

    export TableauRadauIA, TableauRadauIIA

    TableauList = (
       # explicit
       explicit_euler    = TableauExplicitEuler,
       forward_euler     = TableauForwardEuler,
       explicit_midpoint = TableauExplicitMidpoint, 
       heun2             = TableauHeun2,
       heun3             = TableauHeun3,
       ralston2          = TableauRalston2,
       ralston3          = TableauRalston3,
       runge             = TableauRunge,
       runge2            = TableauRunge2,
       kutta             = TableauKutta,
       kutta3            = TableauKutta3,
       rk416             = TableauRK416,
       rk4               = TableauRK4,
       rk438             = TableauRK438,
       ssprk3            = TableauSSPRK3,
       # dirk
       kraaijevanger_spijker = TableauKraaijevangerSpijker,
       qin_zhang             = TableauQinZhang,
       crouzeix              = TableauCrouzeix,
       # firk
       implicit_euler    = TableauImplicitEuler,
       backward_euler    = TableauBackwardEuler,
       implicit_midpoint = TableauImplicitMidpoint,
       crank_nicolson    = TableauCrankNicolson,
       srk3              = TableauSRK3,
    )

    export TableauList

end
