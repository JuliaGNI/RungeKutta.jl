module RungeKutta

    using DelimitedFiles
    using Markdown
    using PrettyTables
    using Reexport
    using StaticArrays

    import GenericLinearAlgebra
    import GeometricBase.Utils: @big, @define
    import LinearAlgebra: istril
    import Polynomials
    import Polynomials: Polynomial


    include("utils.jl")

    include("abstract.jl")
    include("tableau.jl")
    include("tableau_partitioned.jl")

    export Tableau, PartitionedTableau

    export isexplicit,
           isimplicit,
           isdiagnonallyimplicit,
           isfullyimplicit

    export reference

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
    
    export SymplecticTableau,
           SymplecticConjugateTableau,
           SymplecticPartitionedTableau

    export check_symplecticity,
           symplecticity_error,
           issymplectic,
           symplectic_conjugate_coefficients


    include("Tableaus.jl")
    @reexport using .Tableaus

    include("PartitionedTableaus.jl")
    @reexport using .PartitionedTableaus

end
