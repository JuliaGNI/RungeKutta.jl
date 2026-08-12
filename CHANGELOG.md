# Changelog

All notable changes to RungeKutta.jl are documented here. Versions follow
[semantic versioning](https://semver.org) as it applies to Julia's `0.x` series, where a
change to the minor version may break compatibility.

## Unreleased

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
