using RungeKutta
using Test

include("helpers/symtype.jl")

import SymPyPythonCall
using QuadratureRules: gauss_legendre_nodes, gauss_legendre_weights

@testset "$(rpad("Order Conditions",80))" begin
    for s in 1:5
        g = TableauGauss(s)

        for σ in 1:s
            @test check_order_conditions_b(g, σ)
            @test check_order_conditions_c(g, σ) == Array{Bool}(ones(s))
            @test check_order_conditions_d(g, σ) == Array{Bool}(ones(s))
        end

        @test satisfies_simplifying_assumption_b(g)
        @test satisfies_simplifying_assumption_c(g)
        @test satisfies_simplifying_assumption_d(g)
    end

    for s in 2:5
        @test satisfies_simplifying_assumption_b(TableauLobattoIIIA(s))
        @test satisfies_simplifying_assumption_b(TableauLobattoIIIB(s))
        @test satisfies_simplifying_assumption_b(TableauLobattoIIIC(s))
        @test satisfies_simplifying_assumption_b(TableauLobattoIIIC̄(s))
        @test satisfies_simplifying_assumption_b(TableauLobattoIIID(s))
        @test satisfies_simplifying_assumption_b(TableauLobattoIIIE(s))
        @test satisfies_simplifying_assumption_b(TableauRadauIA(s))
        @test satisfies_simplifying_assumption_b(TableauRadauIIA(s))

        @test satisfies_simplifying_assumption_c(TableauLobattoIIIA(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIB(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIC(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIC̄(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIID(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIE(s))
        @test !satisfies_simplifying_assumption_c(TableauRadauIA(s))
        @test satisfies_simplifying_assumption_c(TableauRadauIIA(s))

        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIA(s))
        @test satisfies_simplifying_assumption_d(TableauLobattoIIIB(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIC(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIC̄(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIID(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIE(s))
        @test satisfies_simplifying_assumption_d(TableauRadauIA(s))
        @test !satisfies_simplifying_assumption_d(TableauRadauIIA(s))
    end
end

@testset "$(rpad("Order conditions in BigFloat, Float32 and Float64",80))" begin
    # B(s) holds in exact arithmetic; the target 1/k is formed in the element type, so the
    # residuals lie at the rounding level of BigFloat, below the default tolerance 16eps(BigFloat)
    for f in (TableauGauss, TableauRadauIIA, TableauLobattoIIIA)
        for s in 3:5
            tab = f(BigFloat, s)
            @test eltype(tab) == BigFloat
            @test satisfies_simplifying_assumption_b(tab)
        end
    end

    @test satisfies_simplifying_assumption_b(TableauGauss(Float32, 3))
    @test satisfies_simplifying_assumption_b(TableauGauss(Float64, 3))

    # a caller's explicit tolerance is used as given
    near64 = Tableau(:near64, 1, [0.0;;], [1 + 1e-10], [0.0])
    @test !check_order_conditions_b(near64, 1)
    @test check_order_conditions_b(near64, 1; atol = 1e-9)
    @test check_order_conditions_b(near64, 1; rtol = 1e-9)
end

@testset "$(rpad("Order conditions in Rational{BigInt} and Int",80))" begin
    QQ = Rational{BigInt}

    # the default tolerance of an exact element type is exactly zero, so every explicit
    # tabulated tableau with rational coefficients satisfies B(p) for its order p exactly
    for f in TableauList
        tab = f(QQ)
        isexplicit(tab) || continue
        @test eltype(tab) == QQ
        @test satisfies_simplifying_assumption_b(tab, order(tab))
    end

    rk4 = TableauRK4(QQ)
    @test rk4.s == 4
    @test order(rk4) == 4
    @test check_order_conditions_c(rk4, 1) == trues(4)
    @test check_order_conditions_d(rk4, 1) == trues(4)

    # a weight 1e-30 away from 1 fails B(1) at the exact default tolerance
    near = Tableau(:near, 1, [0//1;;], [1 + 1//big(10)^30], [0//1])
    @test eltype(near) == QQ
    @test !check_order_conditions_b(near, 1)
    @test check_order_conditions_b(near, 1; atol = 1//big(10)^29)
    @test check_order_conditions_b(near, 1; rtol = 1//big(10)^29)

    # an integer element type compares exactly as well
    euler = Tableau(:int, 1, [0;;], [1], [0])
    @test eltype(euler) == Int
    @test check_order_conditions_b(euler, 1)
    @test check_order_conditions_c(euler, 1) == trues(1)
    @test check_order_conditions_d(euler, 1) == falses(1)   # b₁a₁₁ = 0, b₁(1 - c₁) = 1
    @test !check_order_conditions_b(Tableau(:int2, 1, [0;;], [2], [0]), 1)

    # N + 1 and N are equal in Float64, so an integer tableau compares exactly only where
    # the right-hand sides of C and D stay exact and are not converted to Float64
    N = Int64(2)^53
    @test Float64(N + 1) == Float64(N)
    # C(1) reads a₁₁ = c₁
    @test check_order_conditions_c(Tableau(:bigc, 1, [N+1;;], [1], [N]), 1) == falses(1)
    @test check_order_conditions_c(Tableau(:bigc, 1, [N;;], [1], [N]), 1) == trues(1)
    # D(1) reads b₁a₁₁ = b₁(1 - c₁), here with b₁ = 1 and 1 - c₁ = N
    @test check_order_conditions_d(Tableau(:bigd, 1, [N+1;;], [1], [1-N]), 1) == falses(1)
    @test check_order_conditions_d(Tableau(:bigd, 1, [N;;], [1], [1-N]), 1) == trues(1)

    # at the edge of Int64, the difference of the two sides overflows in Int64 rational
    # arithmetic; the checks return false with the default and with explicit tolerances
    M = typemax(Int64)
    for kw in ((;), (; atol = 0, rtol = 0), (; atol = 1))
        # B(1) reads b₁ = 1, here with b₁ = typemin(Int64)
        @test !check_order_conditions_b(Tableau(:edgeb, 1, [0;;], [typemin(Int64)], [0]), 1; kw...)
        # C(1) reads a₁₁ = c₁, here with a₁₁ = M and c₁ = -1
        @test check_order_conditions_c(Tableau(:edgec, 1, [M;;], [1], [-1]), 1; kw...) ==
              falses(1)
        # D(1) reads b₁a₁₁ = b₁(1 - c₁), here with b₁ = 1, a₁₁ = M and 1 - c₁ = -1
        @test check_order_conditions_d(Tableau(:edged, 1, [M;;], [1], [2]), 1; kw...) ==
              falses(1)
    end
end

const QQ = Rational{BigInt}

@testset "$(rpad("Residuals of the simplifying assumptions",80))" begin
    # explicit midpoint rule: a = [0 0; 1/2 0], b = [0, 1], c = [0, 1/2]
    mid = TableauExplicitMidpoint(QQ)
    @test RungeKutta.residual_b(mid, 1) == 0
    @test RungeKutta.residual_b(mid, 2) == 0
    @test RungeKutta.residual_b(mid, 3) == 1 // 4 - 1 // 3
    @test RungeKutta.residual_b(mid, 3) isa QQ
    @test RungeKutta.residuals_c(mid, 1) == [0, 0]
    @test RungeKutta.residuals_c(mid, 2) == [0, -1 // 8]
    @test eltype(RungeKutta.residuals_c(mid, 2)) == QQ
    @test RungeKutta.residuals_d(mid, 1) == [1 // 2, -1 // 2]
    @test eltype(RungeKutta.residuals_d(mid, 1)) == QQ

    for T in (Float64, BigFloat)
        g = TableauGauss(T, 3)
        @test RungeKutta.residual_b(g, 6) isa T
        @test abs(RungeKutta.residual_b(g, 6)) <= 16eps(T)
        @test eltype(RungeKutta.residuals_c(g, 3)) == T
        @test eltype(RungeKutta.residuals_d(g, 3)) == T
    end
end

@testset "$(rpad("order_conditions",80))" begin
    # one entry per tree of order 1 … p, by increasing order, with the residual bᵀψ(t) - 1/γ(t)
    rk4 = TableauRK4(QQ)
    oc = order_conditions(rk4, 5)
    @test length(oc) == 1 + 1 + 2 + 4 + 9
    @test issorted([RungeKutta.tree_order(e.tree) for e in oc])
    @test all(iszero(e.residual) for e in oc if RungeKutta.tree_order(e.tree) <= 4)
    @test any(!iszero(e.residual) for e in oc if RungeKutta.tree_order(e.tree) == 5)
    @test [e.residual for e in oc] ==
          [e.residual for e in order_conditions(rk4, 6)][1:length(oc)]

    # the residuals of the explicit midpoint rule at order 3: bᵀc² - 1/3 and bᵀAc - 1/6
    mid = TableauExplicitMidpoint(QQ)
    oc3 = [e for e in order_conditions(mid, 3) if RungeKutta.tree_order(e.tree) == 3]
    @test sort([e.residual for e in oc3]) == sort([1 // 4 - 1 // 3, 0 - 1 // 6])

    # element types
    for T in (Float32, Float64, BigFloat)
        oc = order_conditions(TableauGauss(T, 2), 4)
        @test eltype(first(oc).residual) == T
        @test all(e -> e.residual isa T, oc)
        @test all(e -> abs(e.residual) <= 16eps(T), oc)
    end
    @test all(e -> e.residual isa QQ, order_conditions(rk4, 4))
    inttab = Tableau{Int}(:int, 1, [0;;], [1], [0])
    @test all(e -> e.residual isa QQ, order_conditions(inttab, 4))
    @test order_conditions(inttab, 1)[1].residual == 0
    @test order_conditions(inttab, 2)[2].residual == -1 // 2

    # SymPy: TableauGauss(symtype(), 2) converts the BigFloat coefficients to Sym; the exact
    # tableau is built from the coefficient functions, as the docs pages do
    @test all(e -> e.residual isa symtype(), order_conditions(TableauGauss(symtype(), 2), 4))
    S = symtype()
    gsym = Tableau(:Gauss, 4, RungeKutta.Tableaus.gauss_coefficients(S, 2),
        gauss_legendre_weights(S, 2), gauss_legendre_nodes(S, 2))
    @test eltype(gsym) == S
    # the residuals of the 2-stage Gauss method vanish to order 4, not at order 5
    ocsym = order_conditions(gsym, 5)
    @test all(e -> e.residual isa symtype(), ocsym)
    @test count(e -> RungeKutta.tree_order(e.tree) <= 4, ocsym) == 8
    @test all(iszero(SymPyPythonCall.simplify(e.residual))
    for e in ocsym if RungeKutta.tree_order(e.tree) <= 4)
    @test any(!iszero(SymPyPythonCall.simplify(e.residual))
    for e in ocsym if RungeKutta.tree_order(e.tree) == 5)

    # p = 0 and p < 0
    @test isempty(order_conditions(rk4, 0))
    @test_throws r"^ArgumentError: the order p must be nonnegative$" order_conditions(rk4, -1)

    # the trees of an entry are copies: changing one leaves the next call unchanged
    e = order_conditions(rk4, 3)[end]
    e.tree.levels[1] = 7
    @test RungeKutta.levels(order_conditions(rk4, 3)[end].tree)[1] == 1
end

@testset "$(rpad("computed_order: tolerance",80))" begin
    # the threshold bracket: the residual of the order-2 tree bᵀc - 1/2 is 3e-10, and its
    # bound rtol · |t| · |b|ᵀ|ψ|(t) = rtol · 2 · (1/2 + 3e-10) is 1.25 and 0.8 times the residual
    tab = Tableau{Float64}(:mid, 2, [0 0; 0.5+3e-10 0], [0, 1], [0, 0.5])
    @test computed_order(tab; rtol = 3.75e-10).order == 2
    @test computed_order(tab; rtol = 2.4e-10).order == 1
    @test linear_order(tab; rtol = 3.75e-10) == 2
    @test linear_order(tab; rtol = 2.4e-10) == 1

    # a small residual that is not round-off
    rk4 = TableauRK4(Float64)
    b = Vector(rk4.b)
    b[1] += 3e-9
    b[2] -= 3e-9
    @test computed_order(Tableau{Float64}(:rk4p, 4, rk4.a, b, rk4.c)).order == 1
    @test computed_order(rk4) == (order = 4, route = :trees, exact = true)

    # an exact default
    rk4q = TableauRK4(QQ)
    bq = Vector(rk4q.b)
    bq[1] += 1 // big(10)^30
    bq[2] -= 1 // big(10)^30
    @test computed_order(Tableau{QQ}(:rk4q, 4, rk4q.a, bq, rk4q.c)).order == 1

    # the default of the other precision, and the cap 2s on B
    @test computed_order(TableauGauss(Float32, 4)).order == 8
    @test computed_order(TableauGauss(Float32, 5)).order == 10
end

@testset "$(rpad("computed_order: route",80))" begin
    @test computed_order(TableauGauss(BigFloat, 3); pmax = 4) ==
          (order = 6, route = :simplifying_assumptions, exact = true)

    # the route reads A𝟙, not the stored c
    a = [0 0 0; 5//24 1//3 -1//24; 1//6 2//3 1//6]
    b = [1 // 6, 2 // 3, 1 // 6]
    for c in ([0, 1 // 2, 1], [0, 1 // 3, 1])
        tab = Tableau{QQ}(:IIIA3, 4, a, b, c)
        @test computed_order(tab; pmax = 2) ==
              (order = 4, route = :simplifying_assumptions, exact = true)
    end

    # the trees of order pmax are tested: an order of pmax - 1 is found by the trees
    @test computed_order(TableauRK4(QQ); pmax = 5) ==
          (order = 4, route = :trees, exact = true)
    @test computed_order(TableauGauss(BigFloat, 2); pmax = 5) ==
          (order = 4, route = :trees, exact = true)

    # the order is a lower bound where the simplifying assumptions do not reach B
    @test computed_order(TableauRK4(QQ); pmax = 3) ==
          (order = 3, route = :trees, exact = false)

    # the bound of the simplifying assumptions: min(q_B, η + ζ + 1, 2η + 2)
    @test RungeKutta._order_bound(4, 3, 3) == 4
    @test RungeKutta._order_bound(8, 3, 1) == 5
    @test RungeKutta._order_bound(8, 1, 5) == 4

    @test_throws r"^ArgumentError: pmax must be at least 1$" computed_order(TableauRK4(QQ);
        pmax = 0)
    # RK4 has B(4), C(1) and D(1): p_BCD = min(4, 3, 4) = 3, above pmax = 1 and below q_B = 4
    @test computed_order(TableauRK4(QQ); pmax = 1) ==
          (order = 3, route = :simplifying_assumptions, exact = false)
end

@testset "$(rpad("computed_order: families",80))" begin
    for (f, srange) in ((TableauGauss, 1:5), (TableauRadauIIA, 2:4), (
        TableauLobattoIIIA, 2:5),
        (TableauLobattoIIIF, 2:5))
        for s in srange
            tab = f(BigFloat, s)
            r = computed_order(tab)
            @test r.order == order(tab)
            @test r.order isa Int
            @test r.exact
            @test r.order <= linear_order(tab)
        end
    end
    @test computed_order(TableauLobattoIIIB(BigFloat, 2)).order == 2
    @test computed_order(TableauLobattoIIIB(Float64, 2)).order == 2
    @test computed_order(TableauExplicitEuler(QQ)) ==
          (order = 1, route = :trees, exact = true)

    for f in TableauList
        tab = f(QQ)
        isexplicit(tab) || continue
        @test computed_order(tab).order == order(tab)
    end

    for s in 2:4, f in (TableauGauss, TableauRadauIIA, TableauLobattoIIIA)

        @test computed_order(f(Float32, s)).order == order(f(Float32, s))
        @test computed_order(f(Float64, s)).order == order(f(Float64, s))
    end
end

@testset "$(rpad("stage_order and isconsistent",80))" begin
    for T in (Float64, BigFloat)
        @test stage_order(TableauGauss(T, 3)) == 3
        @test stage_order(TableauRadauIIA(T, 3)) == 3
        @test stage_order(TableauLobattoIIIA(T, 3)) == 3
        @test stage_order(TableauRadauIA(T, 3)) == 2
        @test stage_order(TableauLobattoIIIC(T, 3)) == 2
        @test stage_order(TableauLobattoIIIB(T, 3)) == 1
        @test stage_order(TableauLobattoIIIB(T, 2)) == 0
        @test stage_order(TableauGauss(T, 3)) isa Int
    end
    @test stage_order(TableauRK4(QQ)) == 1
    @test stage_order(TableauExplicitEuler(QQ)) == 1
    @test stage_order(TableauRK4(QQ); atol = 1, rtol = 0) == 4

    for f in (TableauLobattoIIIB, TableauLobattoIIIĀ, TableauLobattoIIIE, TableauLobattoIIIĒ,
        TableauLobattoIIIF̄, TableauLobattoIIIG, TableauLobattoIIIḠ)
        @test !isconsistent(f(BigFloat, 2))
        @test isconsistent(f(BigFloat, 3))
        @test isconsistent(f(BigFloat, 3)) isa Bool
    end
    for f in TableauList
        tab = f(QQ)
        isexplicit(tab) || continue
        @test isconsistent(tab)
    end
    # c ≠ A𝟙 with the same sum: the first row is 1/2 too large, the second 1/2 too small
    @test !isconsistent(Tableau{QQ}(:swap, 1, [0 0; 1//2 0], [0, 1], [1 // 2, 0]))
    @test isconsistent(Tableau{QQ}(:swap, 1, [0 0; 1//2 0], [0, 1], [1 // 2, 0]); atol = 1)
end

@testset "$(rpad("error_constant",80))" begin
    ec = error_constant(TableauExplicitMidpoint(QQ))
    @test ec isa BigFloat
    @test isapprox(ec, sqrt(big(17)) / 24; rtol = 8eps(BigFloat))
    ec = error_constant(TableauHeun2(QQ))
    @test ec isa BigFloat
    @test isapprox(ec, sqrt(big(5)) / 12; rtol = 8eps(BigFloat))
    @test error_constant(TableauExplicitMidpoint(Float64)) isa Float64
    @test isapprox(error_constant(TableauExplicitMidpoint(Float64)), sqrt(17) / 24;
        rtol = 8eps(Float64))

    @test_throws r"^ArgumentError: the computed order is a lower bound; increase pmax$" error_constant(
        TableauRK4(QQ); pmax = 3)
    @test_throws r"^ArgumentError: pmax must be at least 1$" error_constant(TableauRK4(QQ);
        pmax = 0)

    # the order comes from the simplifying assumptions; the trees of order p + 1 > pmax + 1
    @test computed_order(TableauGauss(QQ, 1); pmax = 1).order == 2
    @test isapprox(error_constant(TableauGauss(BigFloat, 1); pmax = 1),
        error_constant(TableauGauss(BigFloat, 1)); rtol = 8eps(BigFloat))
end

@testset "$(rpad("Order diagnostics: non-finite coefficients and SymPy",80))" begin
    for tab in (Tableau(:bad, 1, [NaN;;], [1.0], [0.0]), Tableau(
        :bad, 1, [0.0;;], [Inf], [0.0]))
        msg = r"^ArgumentError: the tableau has a coefficient that is not finite$"
        @test_throws msg order_conditions(tab, 2)
        @test_throws msg computed_order(tab)
        @test_throws msg error_constant(tab)
        @test_throws msg stage_order(tab)
        @test_throws msg isconsistent(tab)
    end
    # a non-finite node gives no stage order and no consistency, and no error
    @test stage_order(Tableau(:badc, 1, [0.0;;], [1.0], [NaN])) == 0
    @test !isconsistent(Tableau(:badc, 1, [0.0;;], [1.0], [Inf]))

    # the computed order and the error constant need a real element type
    S = symtype()
    gsym = Tableau(:Gauss, 4, RungeKutta.Tableaus.gauss_coefficients(S, 2),
        gauss_legendre_weights(S, 2), gauss_legendre_nodes(S, 2))
    for tab in (gsym, TableauGauss(S, 2))
        @test_throws r"^MethodError: no method matching computed_order\(::(RungeKutta\.)?Tableau\{" computed_order(tab)
        @test_throws r"^MethodError: no method matching error_constant\(::(RungeKutta\.)?Tableau\{" error_constant(tab)
    end
    # SymPy defines eps(Sym), so the checks of C, and with them the stage order, run on Sym
    @test stage_order(gsym) == 2
    @test isconsistent(gsym)
end

@testset "$(rpad("Order diagnostics: degenerate tableaus",80))" begin
    zt = Tableau{QQ}(:zero, 0, zeros(QQ, 2, 2), zeros(QQ, 2), zeros(QQ, 2))
    @test [e.residual for e in order_conditions(zt, 1)] == [-1]
    @test computed_order(zt) == (order = 0, route = :trees, exact = true)
    @test stage_order(zt) == 2
    @test isconsistent(zt)
    @test error_constant(zt) == 1
    @test error_constant(zt) isa BigFloat

    @test computed_order(Tableau{QQ}(:half, 1, [0//1;;], [1 // 2], [0//1])) ==
          (order = 0, route = :trees, exact = true)
end

@testset "$(rpad("Order diagnostics: exports",80))" begin
    for name in (:order_conditions, :computed_order, :stage_order, :error_constant, :isconsistent)
        @test name in names(RungeKutta)
    end
    for name in (:RootedTree, :rooted_trees, :density, :symmetry, :stage_vectors, :residual_b,
        :residuals_c, :residuals_d)
        @test !(name in names(RungeKutta))
    end
end
