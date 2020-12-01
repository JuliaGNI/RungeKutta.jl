module RungeKutta

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

end
