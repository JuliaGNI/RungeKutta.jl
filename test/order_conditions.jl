using RungeKutta
using Test

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
