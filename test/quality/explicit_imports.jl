using ExplicitImports: test_explicit_imports
using RungeKutta

test_explicit_imports(
    RungeKutta;
    # `@big` is imported in RungeKutta for the submodule Tableaus, which reaches it through
    # its parent (`using ..RungeKutta: @big`); ExplicitImports calls it stale (KNOWN_ISSUES.md, K2).
    ignore = (Symbol("@big"),)
)
