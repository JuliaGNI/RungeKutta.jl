module RungeKutta

    include("utils.jl")

    include("tableau.jl")

    export Tableau

    include("order_conditions.jl")

    export check_order_conditions_b,
           check_order_conditions_c,
           check_order_conditions_d,
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

    export TableauCrouzeix

    include("tableaus/erk.jl")

    export TableauExplicitEuler,
           TableauExplicitMidpoint, 
           TableauHeun,
           TableauRunge,
           TableauKutta,
           TableauERK416,
           TableauERK438

    include("tableaus/firk.jl")

    export TableauImplicitEuler,
           TableauImplicitMidpoint,
           TableauSRK3

    include("tableaus/gauss.jl")
    
    export TableauGauss

    include("tableaus/lobatto.jl")

    export TableauLobattoIIIA, TableauLobattoIIIB, TableauLobattoIIIC, TableauLobattoIIIC̄,
           TableauLobattoIIID, TableauLobattoIIIE, TableauLobattoIIIF, TableauLobattoIIIG

    include("tableaus/radau.jl")

    export TableauRadauIIA

end
