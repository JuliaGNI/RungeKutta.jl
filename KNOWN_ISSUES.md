# Known issues

What is known to be broken or incomplete and is not fixed yet. Delete an entry when its issue is
fixed; the fix goes in `CHANGELOG.md`.

### K1 · RungeKutta relies on implicit imports

- **location:** `src/RungeKutta.jl:3`
- **evidence:** `print_explicit_imports(RungeKutta; report_non_public = true)` reports 15 names
  used through an implicit `using`:
  `using DelimitedFiles: DelimitedFiles, readdlm, writedlm`,
  `using Markdown: Markdown`,
  `using PrettyTables: PrettyTables, LatexCell, LatexTableFormat, TextTableBorders, TextTableFormat, pretty_table`,
  `using Reexport: Reexport, @reexport`,
  `using StaticArrays: StaticArrays, SMatrix, SVector`.
- **kind:** defect
- **found:** 2026-08-31

### K2 · ExplicitImports reports `@big` as stale, but it is not

- **location:** `src/RungeKutta.jl:12`
- **evidence:** `check_no_stale_explicit_imports(RungeKutta)` reports `@big` as a stale import of
  module `RungeKutta`. `@big` is imported there from `GeometricBase.Utils`, and the submodule
  `RungeKutta.Tableaus` reaches it through its parent
  (`src/Tableaus.jl:11`: `using ..RungeKutta: big, @big`), then uses it in
  `src/tableaus/erk.jl`, `src/tableaus/dirk.jl` and `src/tableaus/firk.jl`. ExplicitImports does
  not follow that path, so it reports the import as unused in `RungeKutta` itself. Removing the
  import breaks `RungeKutta.Tableaus`.
- **kind:** upstream
- **found:** 2026-09-01
