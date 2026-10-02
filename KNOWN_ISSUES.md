# Known issues

What is known to be broken or incomplete and is not fixed yet. Delete an entry when its issue is
fixed; the fix goes in `CHANGELOG.md`.

### K2 · ExplicitImports reports `@big` as stale, but it is not

- **location:** `src/RungeKutta.jl:13`
- **evidence:** `check_no_stale_explicit_imports(RungeKutta)` reports `@big` as a stale import of
  module `RungeKutta`. `@big` is imported there from `GeometricBase.Utils`, and the submodule
  `RungeKutta.Tableaus` reaches it through its parent
  (`src/Tableaus.jl:11`: `using ..RungeKutta: @big`), then uses it in
  `src/tableaus/erk.jl`, `src/tableaus/dirk.jl` and `src/tableaus/firk.jl`. ExplicitImports does
  not follow that path, so it reports the import as unused in `RungeKutta` itself. Removing the
  import breaks `RungeKutta.Tableaus`.
- **kind:** upstream
- **found:** 2026-09-01

### K3 · No test covers two tableaus with different R∞ in `PartitionedTableau`

- **location:** `src/tableau_partitioned.jl:47`
- **evidence:** the mutant that deletes ` && q.R∞ == p.R∞` survives the whole `core` group
  (`mutate.jl <package> src/tableau_partitioned.jl ' && q.R∞ == p.R∞' '' core`). No test builds a
  `PartitionedTableau` from two tableaus whose R∞ are both present and differ. A test to add:
  `@test ismissing(PartitionedTableau(:x, TableauGauss(1), TableauImplicitEuler()).R∞)`.
- **kind:** missing test
- **found:** 2026-09-28

### K4 · The docs build has its own `text/markdown` heading for SymPy tableaus

- **location:** `docs/make.jl:33-37`
- **evidence:** `docs/make.jl` defines its own
  `Base.show(io, ::MIME"text/markdown", tab::Tableau{symtype()})`. It builds the heading line with
  `Markdown.parse("Runge-Kutta Tableau $(tab.name) with $(tab.s) stages and order $(tab.o):")`,
  a copy of the heading in `src/tableau.jl`. The method in `src/tableau.jl` prints the heading as
  text, so a name with `_` or `*` keeps it literally. The docs method turns such a name into
  emphasis. A fix: print the heading in `docs/make.jl` as `src/tableau.jl` does.
- **kind:** docs
- **found:** 2026-10-02

### K5 · No test pins the line break after the `text/markdown` heading

- **location:** `test/tableau.jl:147-150`
- **evidence:** the tests compare only the first line of `repr("text/markdown", tab)`. The mutant
  that changes the end of the heading in `src/tableau.jl` from `:\n"` to `:\n\n"` survives
  `test/tableau.jl` and `test/quality/explicit_imports.jl`. A test to add:
  ```` @test startswith(repr("text/markdown", gauss), "Runge-Kutta Tableau Gauss with 2 stages and order 4:\n```math\n") ````.
- **kind:** missing test
- **found:** 2026-10-02
