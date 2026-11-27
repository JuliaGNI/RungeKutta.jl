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

- **location:** `docs/make.jl:23-27`
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

### K6 · The real interval of `TableauRK5` in a float type is up to 242 ULPs from its exact value

- **location:** `src/stability.jl:463`
- **evidence:** `real_stability_interval(TableauRK5(Float64))` is `5.603972407468875`. The exact
  interval of the same stored `Float64` tableau, converted to `Rational{BigInt}`, rounds to
  `5.60397240746866`: 242 ULPs apart. In `Float32` the two are `5.6040125f0` and `5.6039724f0`,
  84 ULPs apart. The coefficients of `P` and `Q` carry the round-off of the Faddeev–LeVerrier
  recursion, and the endpoint has the condition number 778. The tests bound the `Float32` value
  within 1024 ULPs.
- **kind:** found late
- **found:** 2026-11-27, the decided rule on the zero counts of the floating-point decisions
  (scenario "Float32 decisions"), which computes `P` and `Q` in the element type of the tableau

### K7 · `isastable` refuses an A-stable tableau whose `P` cancels a pole of `Q` in the left half-plane

- **location:** `src/stability.jl:425`
- **evidence:** `Tableau{Rational{BigInt}}(:red, 1, [-1 0; 0 1], [0, 1], [0, 0])` has
  `P(z) = 1 + z` and `Q(z) = 1 - z²`, so `R(z) = 1/(1 - z)`, which is A-stable.
  `isastable` returns `false`, because condition (i) sees the zero `z = -1` of `Q`. The two
  conditions are equivalent to A-stability only where `P` and `Q` have no common zero.
- **kind:** found late
- **found:** 2026-11-27, condition (i) of the requirement "A- and L-stability by an algebraic test"

### K8 · In `Float32`, the interval of a pair can differ from its value in `Float64`

- **location:** `src/stability.jl:534`
- **evidence:** `imaginary_stability_interval(PartitionedTableau(f(Float32, 5)))` gives
  `13.720513`, `Inf32`, `5.597275` and `18.435873` for Gauss, Lobatto III, Lobatto IIIC and Radau
  IIA, where `imaginary_stability_interval(f(Float32, 5))` gives `Inf32`, `0.0f0`, `Inf32` and
  `Inf32`. `imaginary_stability_interval(TableauLobattoIIIAIIIB(T, 4))` is `5.4772263f0` in
  `Float32` and `3.112717648164224` in `Float64` and `BigFloat`. The tests assert the `Float32`
  pairs of one tableau with itself for `s ≤ 4` only.
- **kind:** found late
- **found:** 2026-11-27, the decided rule on the zero counts of the floating-point decisions,
  whose tested domain is `s ≤ 4`

### K9 · `linear_order` exceeds `deg P + deg Q` at a scale beyond `1/rtol`

- **location:** `src/stability.jl:190`
- **evidence:** for the θ-method
  `Tableau{T}(:θ, 1, [0 0; (1 - θ) θ], [1 - θ, θ], [0, 1])`, `linear_order` returns `4` with
  `θ = 1e8` in `Float32` and with `θ = 1e15` in `Float64`, where `deg P + deg Q = 2`. For
  `k = 3` the residual is `θ² - 1/6`, and the bound `rtol · k · |b|ᵀ|A|^(k−1)𝟙` grows as
  `rtol · 6 θ³`, which exceeds it once `θ` is beyond `1/rtol`.
- **kind:** found late
- **found:** 2026-11-27, the decided bound of the linear-order residual

### K10 · `stability_region` returns `NaN` where `P(z)` and `Q(z)` overflow

- **location:** `src/stability.jl:609`
- **evidence:**
  `stability_region(TableauGauss(2); xlims = (-1e200, 1e200), ylims = (-1.0, 1.0), n = 3).absR[1, 2]`
  is `NaN`; `|R(z)|` there is close to `1`. Both polynomials overflow to infinity in `Float64`,
  and their quotient is `NaN`.
- **kind:** found late
- **found:** 2026-11-27

### K11 · `linear_order(TableauLobattoIIIF(BigFloat, 6))` is `1`

- **location:** `src/stability.jl:190`
- **evidence:** `linear_order(TableauLobattoIIIF(BigFloat, 6))` returns `1`, where the linear
  order of Lobatto IIIF is `2s = 12`. The tests cover `s ≤ 4`. Whether the 256-bit tableau of
  `src/tableaus/lobatto.jl` or the residual bound causes it is not established.
- **kind:** found late
- **found:** 2026-11-27

### K12 · The interval of a pair throws an unclear `ArgumentError` where its conditions overflow

- **location:** `src/stability.jl:353`
- **evidence:** for the θ-method
  `Tableau{Float32}(:θ, 1, Float32[0 0; 1 - θ θ], Float32[1 - θ, θ], Float32[0, 1])` with
  `θ = 1f10`, `imaginary_stability_interval(PartitionedTableau(tab))` throws
  `ArgumentError: invalid rational: zero(BigInt)//zero(BigInt)`. The coefficient `1f20` of `D`
  makes the coefficient `1f40` of `D²` overflow to `Inf32`, and the conversion of the condition's
  coefficients to `Rational{BigInt}` fails. With `θ = 1f8` the result is `Inf32`. The tests do
  not go beyond `θ = 1e8`.
- **kind:** found late
- **found:** 2026-11-27
