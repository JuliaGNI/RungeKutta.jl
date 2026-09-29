using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "ExplicitImports" include("quality/explicit_imports.jl")
    @safetestset "Utility functions" include("utils.jl")
    @safetestset "Tableau" include("tableau.jl")
    @safetestset "Partitioned tableau" include("tableau_partitioned.jl")
    @safetestset "Order conditions" include("order_conditions.jl")
    @safetestset "Symmetry" include("symmetry.jl")
    @safetestset "Symplecticity" include("symplecticity.jl")
    @safetestset "Gauss tableaus" include("tableaus/gauss.jl")
    @safetestset "Lobatto tableaus" include("tableaus/lobatto.jl")
    @safetestset "Radau tableaus" include("tableaus/radau.jl")
    @safetestset "Explicit and implicit tableaus" include("tableaus/tableaus.jl")
    @safetestset "Partitioned tableaus" include("tableaus/prk.jl")
    @safetestset "Tableau list" include("Tableaus.jl")
end
