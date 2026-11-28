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

### K9 · The Lobatto constructors do not infer a concrete tableau type

- **location:** `src/tableaus/lobatto.jl`
- **evidence:** `Base.return_types(TableauLobattoIIIF, (Type{Float64}, Int))` and
  `Base.return_types(TableauLobattoIIIA, (Type{Float64}, Int))` both return
  `Tableau{Float64, _A, Int64} where _A`: the type of the coefficient matrix is not inferred.
- **kind:** defect
- **found:** 2026-11-27

### K10 · The symmetry, symplecticity and order checks overflow for a `Rational{Int64}` tableau near the Int64 limit

- **location:** `src/symmetry.jl:5`
- **evidence:** `check_symmetry`, `check_symplecticity` (`src/symplecticity.jl:5`, `:18`) and the
  left-hand sides of `check_order_conditions_b`, `_c` and `_d` (`src/order_conditions.jl:20`,
  `:25`, `:30`) form their sums in the element type, and checked `Rational{Int64}` arithmetic
  throws. With `M = typemax(Int64)` and `tab = Tableau(:edge, 1, [M // 1;;], [1 // 1], [0 // 1])`,
  both `check_symmetry(tab)` and `check_symplecticity(tab)` raise
  `OverflowError: 9223372036854775807 + 9223372036854775807 overflowed for type Int64`. With
  `tab2 = Tableau(:edge2, 1, [0//1 0//1; 0//1 0//1], [M // 1, M // 1], [0 // 1, 0 // 1])`,
  `check_order_conditions_b(tab2, 1)` raises the same error. `_divide` in
  `src/order_conditions.jl` forms only the right-hand side of the order checks in BigInt: it
  removes the overflow from the subtraction in `isapprox`, not from a sum on the left-hand side.
  A fix: widen integer and rational entries to BigInt before the sums in all five checks.
- **kind:** defect
- **found:** #57

### K11 · The real interval of `TableauRK5` in a float type is up to 242 ULPs from its exact value

- **location:** `src/stability.jl:483`
- **evidence:** `real_stability_interval(TableauRK5(Float64))` is `5.603972407468875`. The exact
  interval of the same stored `Float64` tableau, converted to `Rational{BigInt}`, rounds to
  `5.60397240746866`: 242 ULPs apart. In `Float32` the two are `5.6040125f0` and `5.6039724f0`,
  84 ULPs apart. The coefficients of `P` and `Q` carry the round-off of the Faddeev–LeVerrier
  recursion, and the endpoint has the condition number 778. The tests bound the `Float32` value
  within 1024 ULPs.
- **kind:** found late
- **found:** 2026-11-27, the decided rule on the zero counts of the floating-point decisions
  (scenario "Float32 decisions"), which computes `P` and `Q` in the element type of the tableau

### K12 · `isastable` refuses an A-stable tableau whose `P` cancels a pole of `Q` in the left half-plane

- **location:** `src/stability.jl:442`
- **evidence:** `Tableau{Rational{BigInt}}(:red, 1, [-1 0; 0 1], [0, 1], [0, 0])` has
  `P(z) = 1 + z` and `Q(z) = 1 - z²`, so `R(z) = 1/(1 - z)`, which is A-stable.
  `isastable` returns `false`, because condition (i) sees the zero `z = -1` of `Q`. The two
  conditions are equivalent to A-stability only where `P` and `Q` have no common zero.
- **kind:** found late
- **found:** 2026-11-27, condition (i) of the requirement "A- and L-stability by an algebraic test"

### K13 · In `Float32`, the interval of a pair can differ from its value in `Float64`

- **location:** `src/stability.jl:555`
- **evidence:** `imaginary_stability_interval(PartitionedTableau(f(Float32, 5)))` gives
  `13.720513`, `Inf32`, `5.597275` and `18.435873` for Gauss, Lobatto III, Lobatto IIIC and Radau
  IIA, where `imaginary_stability_interval(f(Float32, 5))` gives `Inf32`, `0.0f0`, `Inf32` and
  `Inf32`. `imaginary_stability_interval(TableauLobattoIIIAIIIB(T, 4))` is `5.4772263f0` in
  `Float32` and `3.112717648164224` in `Float64` and `BigFloat`. The tests assert the `Float32`
  pairs of one tableau with itself for `s ≤ 4` only.
- **kind:** found late
- **found:** 2026-11-27, the decided rule on the zero counts of the floating-point decisions,
  whose tested domain is `s ≤ 4`

### K14 · `linear_order` exceeds `deg P + deg Q` at a scale beyond `1/rtol`

- **location:** `src/stability.jl:190`
- **evidence:** for the θ-method
  `Tableau{T}(:θ, 1, [0 0; (1 - θ) θ], [1 - θ, θ], [0, 1])`, `linear_order` returns `4` with
  `θ = 1e8` in `Float32` and with `θ = 1e15` in `Float64`, where `deg P + deg Q = 2`. For
  `k = 3` the residual is `θ² - 1/6`, and the bound `rtol · k · |b|ᵀ|A|^(k−1)𝟙` grows as
  `rtol · 6 θ³`, which exceeds it once `θ` is beyond `1/rtol`.
- **kind:** found late
- **found:** 2026-11-27, the decided bound of the linear-order residual

### K15 · `stability_region` returns `NaN` where `P(z)` and `Q(z)` overflow

- **location:** `src/stability.jl:630`
- **evidence:**
  `stability_region(TableauGauss(2); xlims = (-1e200, 1e200), ylims = (-1.0, 1.0), n = 3).absR[1, 2]`
  is `NaN`; `|R(z)|` there is close to `1`. Both polynomials overflow to infinity in `Float64`,
  and their quotient is `NaN`.
- **kind:** found late
- **found:** 2026-11-27

### K16 · `linear_order(TableauLobattoIIIF(BigFloat, 6))` is `1`

- **location:** `src/stability.jl:190`
- **evidence:** `linear_order(TableauLobattoIIIF(BigFloat, 6))` returns `1`, where the linear
  order of Lobatto IIIF is `2s = 12`. The tests cover `s ≤ 4`. Whether the 256-bit tableau of
  `src/tableaus/lobatto.jl` or the residual bound causes it is not established.
- **kind:** found late
- **found:** 2026-11-27

### K17 · The interval of a pair throws an unclear `ArgumentError` where its conditions overflow

- **location:** `src/stability.jl:355`
- **evidence:** for the θ-method
  `Tableau{Float32}(:θ, 1, Float32[0 0; 1 - θ θ], Float32[1 - θ, θ], Float32[0, 1])` with
  `θ = 1f10`, `imaginary_stability_interval(PartitionedTableau(tab))` throws
  `ArgumentError: invalid rational: zero(BigInt)//zero(BigInt)`. The coefficient `1f20` of `D`
  makes the coefficient `1f40` of `D²` overflow to `Inf32`, and the conversion of the condition's
  coefficients to `Rational{BigInt}` fails. With `θ = 1f8` the result is `Inf32`. The tests do
  not go beyond `θ = 1e8`.
- **kind:** found late
- **found:** 2026-11-27

### K18 · The stability functions do not support an integer element type

- **location:** `src/stability.jl:166`
- **evidence:** an integer element type is not one of the decided element types (`Float64`,
  `Float32`, `BigFloat`, `Rational{BigInt}`). For `Tableau{Int}`, `stability_at_infinity` returns a
  `Float64`, not an `Int`: `0.5` for `Tableau{Int}(:half, 1, [2;;], [1], [2])` and `Inf` for
  `Tableau{Int}(:ee, 1, [0;;], [1], [0])`. `imaginary_stability_interval(PartitionedTableau(ee))`
  throws `MethodError: no method matching _pmul(::Vector{Int64}, ::Vector{Float64})`
  (`src/stability.jl:237`). The Faddeev–LeVerrier step `-tr(N) / k` (`src/stability.jl:45`)
  divides an integer trace to a `Float64`, which is converted back on its store into the
  coefficient vector of type `Vector{Int}`.
- **kind:** found late
- **found:** #59

### K19 · `real_stability_interval` throws where the endpoint is far below the Cauchy bound of the roots

- **location:** `src/stability.jl:344`
- **evidence:**
  `real_stability_interval(Tableau{Rational{BigInt}}(:k19, 1, [0;;], [big(2)^10001], [0]))`
  throws `the bisection did not locate a root to the precision of BigFloat`. The endpoint is
  `2^-10000`, and the bisection halves the bracket from `[0, 1 + 2^-10000]`, the Cauchy bound of
  the roots, in at most `4 * precision(F) + 4096` steps, fewer than the about `10000 + 256` that
  reach the endpoint to 256 bits. A bisection whose steps scale with the exponent of the root
  would locate it. The tested scales are `1e±8`.
- **kind:** found late
- **found:** #59

### K20 · The touching test of the intervals misses an excursion below the bound between two roots

- **location:** `src/stability.jl:381`
- **evidence:** for a nonzero `rtol`, a sign change ends an interval only where the polynomial
  falls below the bound of its round-off at the midpoint between that root and the next. An
  excursion below the bound away from the midpoint is not seen. With
  `Tableau{Float64}(:rk, 1, [0 0; 1 0], [8/9, 1/9], [0, 1])`,
  `real_stability_interval(tab; rtol = 1/700)` is `8.999999999999996`; with the default `rtol` it
  is `3.0`.
- **kind:** found late
- **found:** #59, the decided rule on the zero counts of the floating-point decisions

### K21 · The intervals of a single tableau run through a common zero of `P` and `Q` on the axis

- **location:** `src/stability.jl:483`
- **evidence:** `Tableau{Rational{BigInt}}(:c1, 1, [-1;;], [0], [0])` has `P = Q = 1 + z`, whose
  common zero `z = -1` is a singular stage system; `real_stability_interval` returns `Inf`.
  `Tableau{Rational{BigInt}}(:axis, 1, [0 1; -1 0], [0, 0], [0, 0])` has `P = Q = 1 + z²`, with
  the common zeros `±i`; `imaginary_stability_interval` returns `Inf` (`src/stability.jl:501`).
  Whether such a point ends an interval is not decided.
- **kind:** found late
- **found:** #59

### K22 · The interval of a pair runs through a singular stage system of the oscillator

- **location:** `src/stability.jl:555`
- **evidence:** for `PartitionedTableau(tab)` with
  `tab = Tableau{Rational{BigInt}}(:axis, 1, [0 1; -1 0], [0, 0], [0, 0])`, the denominator is
  `D(ν) = (1 - ν²)²`, so the stage system of one step is singular at `ν = 1`, and
  `N1 = D`, `N2 = 4D`. `imaginary_stability_interval(PartitionedTableau(tab))` returns `Inf`, in
  `Rational{BigInt}` and in `Float64`.
- **kind:** found late
- **found:** #59

### K23 · The box scan of the tests skips the callable of `StabilityFunction`

- **location:** `test/stability.jl:618`
- **evidence:** the testset "No captured variable is boxed" keeps the names of `RungeKutta` whose
  value is a `Function`. `StabilityFunction <: Function` is `false`, so the method
  `(R::StabilityFunction)(z)` at `src/stability.jl:24` is not scanned. It holds no closure and no
  box.
- **kind:** missing test
- **found:** #59

### K24 · Some `BigFloat` tableaus are less accurate than `eps(BigFloat)`, so `computed_order` and `stage_order` miss their orders

- **location:** `src/tableaus/radau.jl:61`
- **evidence:** `maximum(abs.(tab.a * ones(BigFloat, s) .- tab.c)) / eps(BigFloat)` is 39.6, 129
  and 262 for `TableauRadauIA(BigFloat, s)`, s = 4, 5, 6, 92.1 and 272 for
  `TableauLobattoIIIB(BigFloat, s)`, s = 5, 6, and 58.5 for `TableauLobattoIIIF(BigFloat, 6)`
  (`src/tableaus/lobatto.jl:265` and `:480`). The largest residual of a tree condition over its
  bound `16eps(BigFloat) · |t| · |b|ᵀ|ψ|(t)` is 1.38 for Radau IA(5), 2.81 for Lobatto IIIF(6) and
  0.985 for Radau IB(5). So `computed_order` gives 5 for Radau IA(5), whose order is 9, 1 for
  Lobatto IIIF(6), whose order is 10, and `(order = 10, route = :trees, exact = false)` for
  Radau IA(6), whose order is 11; `stage_order` gives 0 for Radau IA(4), (5) and (6) and Lobatto
  IIIB(5) and (6). In `Float64` the same tableaus are within 0.5 eps, and
  `computed_order(TableauRadauIA(Float64, 6))` is 11. The tests use `BigFloat` only up to
  Gauss(5), Radau IIA(4), Lobatto IIIA(5) and IIIF(5).
- **kind:** found late
- **found:** 2026-11-28; the same inaccuracy as K16

### K25 · In `Float32`, `computed_order` gives one order more than the true order for some tableaus with six stages

- **location:** `src/order_conditions.jl:217`
- **evidence:** `computed_order(f(Float32, 6))` is
  `(order = 12, route = :simplifying_assumptions, exact = true)` for `TableauRadauIA` and
  `TableauRadauIIA` (order 11), and the same with order 11 for `TableauLobattoIIIA` and
  `TableauLobattoIIIC` (order 10). The count of B on these tableaus is 12 and 11, one more than
  in exact arithmetic: the residual of B(12) and B(11), 3.9e-7 for Radau and 1.7e-6 for Lobatto
  (the same weights and nodes in IIIA and IIIC) in `BigFloat`, lies below the absolute tolerance `16eps(Float32)` = 1.9e-6.
  The counts of C and D are the exact ones. In `Float64` each gives its order. The tests use
  `Float32` only up to four stages.
- **kind:** found late
- **found:** 2026-11-28
