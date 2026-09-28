using ExplicitImports: test_explicit_imports
using RungeKutta

test_explicit_imports(
    RungeKutta;
    # RungeKutta relies on implicit imports from five packages (KNOWN_ISSUES.md, K1).
    no_implicit_imports = false,
    # RungeKutta imports non-public names from GeometricBase.
    all_explicit_imports_are_public = false,
    # RungeKutta qualifies non-public names of GeometricBase and `Core.eval`; on Julia 1.10,
    # which has no `public`, also the names of LinearAlgebra.
    all_qualified_accesses_are_public = false,
    # `@big` is imported in RungeKutta for the submodule Tableaus, which reaches it through
    # its parent (`using ..RungeKutta: @big`); ExplicitImports calls it stale (KNOWN_ISSUES.md, K2).
    ignore = (Symbol("@big"),)
)
