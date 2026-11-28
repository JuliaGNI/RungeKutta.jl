# Changelog

All notable changes to RungeKutta.jl are documented here. Versions follow
[semantic versioning](https://semver.org) as it applies to Julia's `0.x` series, where a
change to the minor version may break compatibility.

## Unreleased

### Added

- The stability function and the linear and nonlinear stability of a tableau, in the new file
  `src/stability.jl`. `stability_function(tab)` returns a `StabilityFunction`, the coefficients of
  `R(z) = P(z)/Q(z)` in the element type of the tableau, exact for `Rational{BigInt}` and SymPy,
  and callable as `R(z)`.
  From it, `stability_at_infinity`, `linear_order` and `pade_type` give `R(∞)`, the linear order
  (the order of `R` as an approximation of `exp`, which can differ from the classical order) and
  the Padé type; `isastable` and `islstable` decide A- and L-stability from the coefficients,
  without sampling; `real_stability_interval` and `imaginary_stability_interval` give the
  stability intervals on the two axes, and the second also the interval of a `PartitionedTableau`
  on the harmonic oscillator; `algebraic_stability_matrix` and `isalgebraicallystable` give
  algebraic stability; and `stability_region` returns the values of `|R|` on a grid, for a plot.
  Each predicate takes a keyword `rtol`, whose default is zero for a rational element type, so
  that the decisions are exact there. `stability_at_infinity` computes `R(∞)` and does not read the
  tableau's stored `R∞`. The stability intervals throw an error where the bisection does not
  locate an endpoint to the precision of `float(T)`, rather than return an inexact endpoint.
- The order conditions and the classical order of a tableau, from the rooted trees of Butcher's
  theory, in the new file `src/rooted_trees.jl` and in `src/order_conditions.jl`.
  `order_conditions(tab, p)` returns the residual `bᵀψ(t) - 1/γ(t)` of each tree of order 1 to
  `p`, with no tolerance, in the element type of the tableau, exact for `Rational{BigInt}` and
  SymPy. `computed_order(tab; pmax = 10)` tests these conditions up to order `pmax`, and above
  it takes the order that the simplifying assumptions B, C and D give; its result names the
  route and whether the order is exact or a lower bound. `stage_order(tab)` gives the largest
  `q ≤ s` with C(1) to C(q), `error_constant(tab)` the 2-norm of the coefficients of the leading term
  of the local error, and `isconsistent(tab)` tells whether the stored nodes are `c = A𝟙`. The
  trees use `A𝟙`, not the stored `c`. The test suite compares the trees, their densities and
  symmetries and the elementary weights with RootedTrees.jl, a test dependency only.
- `stabilityplot(tab)` and `stabilityplot!(ax, tab)` plot the stability region of a tableau: the
  set where `|R(z)| ≤ 1` filled, its boundary, and the real and imaginary axes. They are a Makie
  recipe in the package extension `RungeKuttaMakieExt`, which loads with Makie, for example through
  a backend such as CairoMakie; Makie is a weak dependency. The data come from `stability_region`,
  with the same keywords `xlims`, `ylims` and `n` and the same defaults. The axis limits follow the
  window, also when it excludes 0.

## v0.6.8

### Fixed

- `order` of `TableauLobattoIIIF(s)`, `TableauLobattoIIIF̄(s)`, `TableauLobattoIIIG(s)` and
  `TableauLobattoIIIḠ(s)` returns `2s - 2`, their classical order, where it returned `2s`. The
  `s` Lobatto weights satisfy the simplifying assumption B(2s - 2) and not B(2s - 1), so no
  tableau with these weights reaches order `2s`. The pairs `TableauLobattoIIIFIIIF̄(s)`,
  `TableauLobattoIIIF̄IIIF(s)` and `TableauLobattoIIIGIIIḠ(s)` take the smaller order of their
  tableaus, so their `order` changes from `2s` to `2s - 2` as well. The stability function of
  Lobatto IIIF and IIIF̄ is the (s, s)-Padé approximation of `exp(z)`, so their linear order is
  `2s`, which their docstrings state beside the classical order. The stored `R∞` of every
  tableau is unchanged.
- `check_order_conditions_b` compares with `one(T) / k` in the element type `T` of the tableau,
  where it compared with the `Float64` value `1/k`. `satisfies_simplifying_assumption_b` therefore
  holds for the Gauß, Radau IIA and Lobatto IIIA tableaus with 3, 4 and 5 stages in `BigFloat`,
  where it returned `false`.
- The default `atol` and `rtol` of `check_order_conditions_b`, `check_order_conditions_c` and
  `check_order_conditions_d`, and so of the three `satisfies_simplifying_assumption_*`, are
  exactly zero for an integer or rational element type, and these checks form the right-hand
  side `x / k` as the exact rational `x // k` in `BigInt` for such a type, so they compare
  exactly where the sums of the left-hand side do not overflow the element type. They
  raised a `MethodError` for `Rational{BigInt}`, because `eps` has no method for it. Other element
  types keep the default `16eps(T)` and the division in the element type.
- The default `atol` and `rtol` of `check_symplecticity` (for a `Tableau` and a
  `PartitionedTableau`) and of `check_symmetry` are exactly zero for an integer or rational
  element type, so `issymplectic` and `issymmetric` decide exactly for such a tableau where its
  sums do not overflow the element type, as checked `Rational{Int64}` arithmetic can near the
  Int64 limit. In 0.6.7
  all four raised a `MethodError` for `Rational{BigInt}`, because `eps` has no method for it. For
  `Rational{BigInt}`, `issymplectic` and `issymmetric` are now true for the implicit midpoint rule
  and false for `TableauRK4`. `check_symmetry` and `issymmetric` of a `PartitionedTableau` forward
  to the `Tableau` method and are fixed with it. Other element types keep the default `16eps(T)`.
- The methods page lists `TableauRK4`, an alias of `TableauRK416`, with 4 stages and order 4, where
  it listed 6 stages and order 5, and Lobatto IIIF and IIIG with order `2s-2`. The reference of
  Lobatto IIIF names its journal, the IAENG International Journal of Applied Mathematics,
  volume 46, issue 1, and no longer gives page numbers, which the paper does not print.
- The page on Runge-Kutta methods states the quadrature conditions of the weights as B(2s) for
  the Gauß methods, B(2s-1) for the Radau methods and B(2s-2) for the Lobatto methods, where it
  stated B(s), B(s-1) and B(s-2).

## v0.6.7

### Changed

- The documentation no longer uses Weave. The pages of the Gauß, Radau I, Radau II and Lobatto III
  tableaus are Markdown sources with Documenter `@eval` blocks, which print the same symbolic
  tableaus. This fixes the `Documentation` job, which failed while it loaded SymPyPythonCall
  (known issue K6): Weave 0.10.12 holds JSON at 0.21, so the docs environment could not install
  CondaPkg 0.2.38 and PythonCall 0.9.37, and PythonCall 0.9.36 finds no shared library in
  conda-forge's Python 3.15. No change to the package itself.

### Fixed

- `TableauRK42` has order four (issue #54). Its coefficient `a[4,3]` was `1`, where the order
  conditions require `2`, so the row sums of `a` did not equal `c`, four order conditions of
  order three and four failed, and an integrator with this tableau converged with order one.
  A new test checks the row sums and the order conditions up to order four of every explicit
  tableau, and `scripts/verify_rk42_order.jl` checks the tableau in exact rational arithmetic
  and by its observed convergence order.

## v0.6.6

### Changed

- The `StaticArrays` compat floor rises from `1` to `1.5.17` (the upper bound `< 2` is unchanged).
  The General registry caps the Julia compatibility of StaticArrays 1.0.0 through 1.5.16 below
  Julia 1.12, so on this package's Julia floor no earlier StaticArrays can be installed, and the
  `Downgrade` job's forced-minimum check failed. This is a compatibility change only, with no
  behaviour change: a caller who could install this package on Julia 1.12 already had StaticArrays
  1.5.17 or later.
- The `PrettyTables` compat floor rises from `3` to `3.0.9` (the upper bound `< 4` is
  unchanged). With PrettyTables 3.0.0, `show` of a `Tableau` throws `MethodError: no method
  matching length(::Nothing)`, because the text backend fails on a table without column labels,
  which is how `show` prints the coefficients. PrettyTables 3.0.9 fixes this error (PrettyTables
  issue #270).
- The `Aqua` compat floor of the test environment rises from `0.8` to `0.8.14`. Earlier 0.8
  releases read `Core.TypeName.mt` in `test_piracies`, a field that Julia 1.12 does not have,
  so `Aqua.test_all` errors on this package's Julia floor; 0.8.14 reads `Core.methodtable`
  instead.

## v0.6.5

### Changed

- **Julia 1.12 is the floor** (`[compat] julia = "1.12"`, from 1.11). A patch release, because no
  API changes: a caller on Julia 1.11 can no longer install this version.
  The `text/markdown` heading of `Tableau` keeps printing as text, and Markdown stays out of the
  dependencies, as in v0.6.4.
- CI uploads coverage from the `Julia 1 - ubuntu-latest` job instead of `Julia min`, and a
  test job saves the Julia cache only when it succeeds.

## v0.6.4

### Changed

- The Julia floor returns to 1.11, and the Markdown dependency goes. The `text/markdown`
  `show` method of `Tableau` prints its heading line as text instead of passing it through
  `Markdown.parse`, so no qualified access to `Markdown.parse` is left, and that function is
  public only from Julia 1.12.
- The `text/markdown` heading of a tableau whose name has underscores keeps the name
  literally. Before, `Markdown.parse` turned a name such as `my_tab_name` into emphasis.

## v0.6.3

### Changed

- Test suite reorganised: files now mirror `src/` structure, dependencies moved to
  `test/Project.toml`, per-file `@safetestset` groups introduced with filtering via
  `GROUPS` from `ARGS`, and Aqua checks added. No tests lost.
- Project.toml now carries compat entries `LinearAlgebra = "1"` and `Markdown = "1"`, so
  Aqua's deps_compat check passes (issue #32).
- The unused `big` import from `RungeKutta` in `src/Tableaus.jl` is removed. The module
  only uses the `@big` macro (from GeometricBase.Utils), which expands numeric literals
  at macro-expansion time using its own internal helper; the plain `big` binding was stale.
  Bit-identity of all 22 nullary tableau constructors (Float64 and BigFloat) confirmed by
  `scripts/checksum_tableaus.jl` on the base commit and this branch in separate Julia
  processes; the full test suite (2,543 tests) passes, confirming nothing else was affected.
- `test/order_conditions.jl` asserts `check_order_conditions_b(g, σ)` directly, not
  `== true`, so a non-`Bool` result fails the test. No source file changes.
- `test/Project.toml` no longer carries `QuadratureRules = "0.2.1"` in `[compat]`. A test or
  docs environment gives no `[compat]` entry for a dependency of the root `Project.toml`: the
  resolver already applies the root's bound (`"0.2"`), and a second entry can only narrow it.
- The test suite checks explicit imports with ExplicitImports.jl
  (`test/quality/explicit_imports.jl`, `core` group). `@big` is in its `ignore`, because the
  submodule `Tableaus` reaches it through `RungeKutta`.
- `fatou lint` reports no finding in `src/` or `test/`. The nested `if` in the
  `PartitionedTableau` constructor is one condition, and `from_file` no longer assigns the
  unused variable `T`.
- The floors rise to Julia 1.12, GeometricBase 0.15.0, QuadratureRules 0.2.2 and
  CompactBasisFunctions 0.4.2, because GeometricBase 0.15 declares its stubs public and
  requires Julia 1.11. The Julia floor is 1.12, not 1.11, because `Markdown.parse` is public
  only from Julia 1.12, and the explicit-imports check requires every qualified access to be
  public.
- The five `using` lines of `src/RungeKutta.jl` name every symbol that `src/` takes from
  DelimitedFiles, Markdown, PrettyTables, Reexport and StaticArrays.
  `test/quality/explicit_imports.jl` runs every check of `test_explicit_imports`; only `@big`
  stays in its `ignore`.

### Fixed

- `from_file` no longer evaluates the third word of a tableau file's header in `Main`. It ran
  code from the file and discarded the result. A header whose third word does not parse, or
  names an undefined name, no longer throws.
- The documentation builds again. `docs/make.jl` passed the Weave sources as paths relative to
  the working directory, so it worked only when run from `docs/`. The shared Documentation
  workflow runs it from the repository root, and every docs build since v0.6.1 stopped at
  `opening file "src/gauss.jmd": No such file or directory`. The paths now resolve from
  `docs/make.jl` itself.

## v0.6.2

### Changed

- Every tracked file is now Unicode NFC-normalised. Thirteen stored `ā` (38 times), `Ḡ` (28), `Ē`
  (28), `Ā` (27), `â` (9), `ĉ` (9), `Ã` (6), `Â` (5), `é` (4) and `ã` (3) as a base letter plus a
  combining mark, inherited from macOS rather than chosen. Eleven are `.jl`; the other two are the
  Weave documents `docs/src/radau1.jmd` and `docs/src/radau2.jmd`, which `docs/make.jl` executes.

  Nothing about the compiled code changes: Julia's parser normalises identifiers to NFC, so the
  symbols were already precomposed and dispatch is untouched. What changes is that the source now
  matches what a keyboard, an editor search or a `grep` pattern produces — in an NFD file a pattern
  typed in NFC matches nothing at all, silently.

  **Two kinds of string literal change at runtime**, since string literals are *not*
  parser-normalised.

  `Symbol("LobattoIIIAIIIĀ", s)`, `Symbol("LobattoIIIEIIIĒ", s)` and `Symbol("LobattoIIIGIIIḠ", s)`
  in `src/tableaus/prk.jl` now produce precomposed symbols where they produced decomposed ones.
  These become the `name` field of the returned `PartitionedTableau` — shown, and also compared in
  `Base.isequal` (`src/tableau_partitioned.jl:94`), though not in `==` or `hash`. Both operands
  there are constructed at run time, so both are NFC and the comparison is unaffected. Nothing in
  this package or downstream compares such a name against a symbol written in source, so no result
  changes. Code that did compare, having had to spell the name decomposed to match, would need the
  ordinary spelling now.

  The `reference(::Val{…})` bodies in `src/tableaus/lobatto.jl` and `src/tableaus/radau.jl` also
  change: `Padé` is now precomposed. `reference` is exported, so a caller sees the new bytes. The
  tests compare `reference` against `reference` rather than against a literal
  (`test/test_lobatto.jl:350`, `test/test_radau.jl:108–111`), and `Padé` appears nowhere else in
  the tree.

  Note that the file looked mixed before — `Symbol("LobattoIIIBIIIB̄", s)` and its `C̄`, `D̄` and `F̄`
  siblings were already NFC — but that was not an inconsistency: a macron over `B`, `C`, `D` or `F`
  has no precomposed codepoint, so NFC leaves those decomposed and they are unchanged here.

## v0.6.1

### Fixed

- **The `lobatto_nullvector` docstring no longer breaks downstream documentation builds.**
  It referred to the internal helper as ``[`_nullvector`](@ref)``. `_nullvector` lives in
  `RungeKutta`, not in `RungeKutta.Tableaus`, so the link resolves only in a build that
  documents the parent module too, as this package's own `docs/src/library.md` does. A
  downstream package that pulls the accessors in with `@autodocs Modules =
  [RungeKutta.Tableaus]` inherits the docstring but not the target, and Documenter then
  fails the build with a `:cross_references` error. Reported from GeometricIntegrators.jl,
  whose documentation job this broke on the upgrade to v0.6.0. The reference is now plain
  code formatting; nothing about the function itself changed.

## v0.6.0

This release is **breaking**: the tableau accessors are renamed, one of them is no longer
exported, and the dependency bounds are tightened. Downstream code has to be updated; there
are no deprecated aliases for the old names. In JuliaGNI, GeometricIntegrators.jl is the
package affected.

### Breaking

- **The `get_` prefix is dropped from the twelve tableau accessors that remain.**
  `get_gauss_coefficients` becomes `gauss_coefficients`, `get_lobatto_a_coefficients` becomes
  `lobatto_a_coefficients`, `get_lobatto_nullvector` becomes `lobatto_nullvector`, and so on
  for the Gauss, Lobatto and Radau coefficients and the Lobatto nullvector. The other eight of
  the twenty are removed outright, see below. The prefix carried no information, and the
  unprefixed names read consistently with `nodes`, which is now extended from GeometricBase.
  The old names are **removed**, not deprecated, so callers fail with an `UndefVarError`
  rather than on a warning they might not see.
- **`lobatto_nullvector` is no longer exported.** It was the only accessor that was, which
  made it the odd one out; the remaining twelve are now reached alike through
  `RungeKutta.Tableaus`, either qualified or by explicit import.
- **`CompactBasisFunctions` is restricted to `0.3`** and **`GeometricBase` to `0.14.8`**,
  dropping the previously allowed `0.2` and `0.10`–`0.13` respectively.
- **`QuadratureRules` `0.2` is now a dependency.** It is where the quadrature nodes and
  weights come from, so it is required, not optional.
- **The eight node and weight accessors are removed entirely**, not renamed:
  `get_gauss_nodes`, `get_gauss_weights`, `get_lobatto_nodes`, `get_lobatto_weights`,
  `get_radau_1_nodes`, `get_radau_1_weights`, `get_radau_2_nodes`, `get_radau_2_weights`.
  There is no `gauss_nodes` or `lobatto_weights` to migrate to. Once QuadratureRules learned to
  evaluate its nodes and weights on symbolic element types too, these had become one-line
  forwardings, and two public names for one function only invite drift. Call
  `QuadratureRules.gauss_legendre_nodes(T, s)`, `lobatto_legendre_weights(T, s)`,
  `radau_legendre_nodes(T, s, Val(:left))` for Radau IA and `Val(:right)` for Radau IIA.
  The QuadratureRules names are also the more precise ones, `gauss_nodes` being ambiguous
  now that Gauss-Chebyshev nodes exist alongside Gauss-Legendre. Note that the one-argument
  forms defaulted to `BigFloat` whereas QuadratureRules defaults to `Float64`, so pass the
  element type explicitly. `nodes(tab)`, `weights(tab)` and `coefficients(tab)` remain the
  way to read a tableau that has already been constructed.
- **The `normalize` keyword of `lobatto_nullvector` is removed.** It never controlled
  normalisation — the vector was of unit length either way — only whether the sign was fixed.
  Both are now unconditional, so the keyword had nothing left to select. Callers that passed
  `normalize=true` get the same vector by dropping the argument.
- **A one-stage Radau rule no longer throws.** The nodes and weights of the one-node Radau
  quadrature are perfectly well defined — it is a Riemann sum — and QuadratureRules returns
  them. What remains undefined is the one-stage Radau *tableau*, so
  `radau_1_coefficients(1)`, `radau_2_coefficients(1)` and the four `TableauRadau*(1)`
  constructors still throw an `ErrorException`.

### Added

- `coefficients`, `nodes` and `weights` are extended from GeometricBase and exported, so a
  tableau's entries are reached through the same generic functions as elsewhere in JuliaGNI.
  Previously all three were defined bare, which created RungeKutta-local functions: with any
  other package of the ecosystem in scope, `weights(tab)` resolved to that package's generic
  and raised a `MethodError`.

### Changed

- **The Gauss, Lobatto and Radau nodes and weights are taken from QuadratureRules.jl**
  instead of being computed here. QuadratureRules works on `[-1,+1]` and shifts afterwards,
  where the monomial basis is far better conditioned than for the polynomials shifted to
  `[0,1]` that were used before, so the values are more accurate. The worst violation of the
  moment conditions at `s = 10`, measured in `BigFloat`, falls from `5.7e-68` to `9.8e-75`
  for Gauss and from `3.6e-73` to `4.0e-77` for Lobatto; the residual of the defining
  polynomial at the Radau nodes improves by up to a factor 800. Results therefore differ from
  previous releases in the last digits.
- The Radau weights are now evaluated from the Radau closed form rather than by solving the
  Vandermonde system of the simplifying assumption `B(s)`, which is better conditioned for
  many stages.
- The delegation is unconditional in the element type. QuadratureRules evaluates its nodes and
  weights on symbolic element types too, so the tableaux still come out in closed form for
  those, as the manual shows them, without this package keeping a fallback of its own.
- `lobatto_nullvector` returns a vector of unit length whose first entry is positive. The
  nullvector is now spanned by a column-pivoted QR factorisation instead of an SVD, and its
  sign would otherwise have followed the pivot order: for `s = 2` and `s = 5` it came out
  opposite to the sign of previous releases.

### Removed

- **The `GenericLinearAlgebra` dependency.** It was never called by name; it was loaded for
  its methods, which supplied the two things the standard library does not provide for
  `BigFloat`: the eigenvalues of the companion matrix behind `Polynomials.roots`, and the
  singular value decomposition behind `nullspace`. Neither is needed any more — the nodes
  come from QuadratureRules, and the Lobatto nullvector is now obtained from a
  column-pivoted QR factorisation, which is implemented generically in `LinearAlgebra`.
- **The `Polynomials` dependency.** With the nodes and weights coming from QuadratureRules for
  every element type, the last callers of `_legendre` and `_shifted_legendre` were gone, and
  those two helpers were all this package used it for. Both are removed with it.

### Fixed

- **RungeKutta.jl again precompiles on Julia 1.13.** Since v0.4.0, GenericLinearAlgebra
  defines `LinearAlgebra.eigencopy_oftype` for `UpperHessenberg` under a `VERSION < v"1.14"`
  guard. Julia 1.13 added that method to LinearAlgebra itself, so loading the package there
  overwrote an existing method, which is an error during precompilation. Raising the Julia
  bound would not have helped: 1.13 added only the `UpperHessenberg` wrapper plumbing, and
  every dense `eigvals!` and `svd!` remains restricted to the BLAS element types.

## v0.5.23

### Added

- `TableauIRK3`, the two-stage, 3rd order fully implicit Runge-Kutta tableau.

## v0.5.22

### Fixed

- The order of `TableauKraaijevangerSpijker` is corrected from 2 to 1. The tableau itself is
  unchanged; only the reported order was wrong. Code that dispatches on or asserts against
  `order(TableauKraaijevangerSpijker())` sees a different value than in earlier releases.

## v0.5.21

### Changed

- The symbolic backend used by the tests and the documentation moves from `SymPy` to
  `SymPyPythonCall`, replacing the PyCall bridge with PythonCall and CondaPkg. This affects
  only the test and documentation environments, not the package API: `src/` contains no
  reference to SymPy and is generic in the coefficient element type. The practical difference
  is that the Python environment is now provisioned reproducibly by CondaPkg instead of
  requiring a manually configured system Python.

## v0.5.16 – v0.5.20

### Changed

- Dependency bounds and CI housekeeping: compat updates from CompatHelper, a Julia compat
  floor of 1.10, and CI coverage extended through Julia 1.12.
