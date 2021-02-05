module RungeKutta

    using DelimitedFiles
    using Markdown
    using PrettyTables
    using Reexport

    import GenericLinearAlgebra
    import LinearAlgebra: istril
    import Polynomials
    import Polynomials: Polynomial


    include("utils.jl")

    include("tableau.jl")
    include("tableau_partitioned.jl")

    export Tableau, PartitionedTableau

    export isexplicit,
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


    include("Tableaus.jl")
    @reexport using .Tableaus

    include("PartitionedTableaus.jl")
    @reexport using .PartitionedTableaus

end
