using RungeKutta
using Test

import LinearAlgebra: I
import SymPyPythonCall

include("helpers/symtype.jl")

const RB = Rational{BigInt}

# The tableaus that the tests build by hand, each in the element type `T` of its scenario.
θmethod(::Type{T}, θ) where {T} = Tableau{T}(:θ, 1, [0 0; (1 - θ) θ], [1 - θ, θ], [0, 1])
zerotableau(::Type{T}) where {T} = Tableau{T}(:zero, 1, [0;;], [0], [0])
poletableau(::Type{T}) where {T} = Tableau{T}(:pole, 1, [-1;;], [-1], [-1])
inconsistent(::Type{T}) where {T} = Tableau{T}(:inconsistent, 1, [0;;], [1 // 2], [0])
function rkc2(::Type{T}) where {T}
    Tableau{T}(:rkc2, 1, [0 0; 1//4 0], [1 // 2, 1 // 2], [0, 1 // 4])
end

# The coefficients of the (k, j)-Padé approximation of exp(z), padded with zeros to length n.
function pade(k, j, n)
    f = factorial ∘ big
    p = [f(k + j - i) * f(k) // (f(k + j) * f(i) * f(k - i)) for i in 0:k]
    q = [(-1)^i * f(k + j - i) * f(j) // (f(k + j) * f(i) * f(j - i)) for i in 0:j]
    return vcat(p, zeros(RB, n - k - 1)), vcat(q, zeros(RB, n - j - 1))
end

# R(z) = 1 + z bᵀ(I − zA)⁻¹𝟙, solved in the complex type of z.
function resolvent(tab, z)
    A = Matrix{typeof(z)}(tab.a)
    b = Vector{typeof(z)}(tab.b)
    return 1 + z * sum(b .* ((I - z * A) \ ones(typeof(z), tab.s)))
end

# The root of the polynomial with ascending coefficients c, by Newton's method at 512 bits.
function newton_root(c, x0)
    setprecision(BigFloat, 512) do
        p = BigFloat.(c)
        dp = [k * p[k + 1] for k in 1:(length(p) - 1)]
        x = BigFloat(x0)
        for _ in 1:200
            x -= evalpoly(x, p) / evalpoly(x, dp)
        end
        x
    end
end

# |x - y| within n ULPs, in the coarser of the precisions of x and y (a 512-bit reference keeps its
# precision under `BigFloat(y)`).
within(x, y, n) = abs(x - y) <= n * max(eps(abs(x)), eps(abs(y)))

# x has the type and the value of v; `===` compares two `BigFloat`s by identity, not by value.
isvalue(x, v) = x isa typeof(v) && x == v

const families_lo = (TableauRadauIA, TableauRadauIIA, TableauRadauIB, TableauLobattoIII,
    TableauLobattoIIIA, TableauLobattoIIIC, TableauLobattoIIIF, TableauLobattoIIIG)

@testset "$(rpad("The stability function",80))" begin
    # closed forms in exact arithmetic
    for (tab, P, Q) in (
        (TableauExplicitEuler(RB), [1, 1], [1, 0]),
        (TableauImplicitEuler(RB), [1, 0], [1, -1]),
        (TableauImplicitMidpoint(RB), [1, 1 // 2], [1, -1 // 2]),
        (TableauCrankNicolson(RB), [1, 1 // 2, 0], [1, -1 // 2, 0]),
        (TableauRK4(RB), [1, 1, 1 // 2, 1 // 6, 1 // 24], [1, 0, 0, 0, 0]))
        R = stability_function(tab)
        @test R isa StabilityFunction{RB}
        @test R.P isa Vector{RB}
        @test R.Q isa Vector{RB}
        @test length(R.P) == length(R.Q) == tab.s + 1
        @test R.P == P
        @test R.Q == Q
    end

    # the Padé coefficients of the collocation families
    for s in 1:4,
        (family, k, j) in ((TableauGauss, s, s), (TableauRadauIIA, s - 1, s),
            (TableauLobattoIIIA, s - 1, s - 1), (TableauLobattoIIIC, s - 2, s))

        family == TableauGauss || s ≥ 2 || continue
        R = stability_function(family(BigFloat, s))
        P, Q = pade(k, j, s + 1)
        @test R isa StabilityFunction{BigFloat}
        @test all(abs.(R.P .- P) .<= 1e-70)
        @test all(abs.(R.Q .- Q) .<= 1e-70)
    end
    for T in (Float64, Float32), s in 1:4

        R = stability_function(TableauGauss(T, s))
        P, Q = pade(s, s, s + 1)
        @test R isa StabilityFunction{T}
        @test all(within.(R.P, T.(P), 16))
        @test all(within.(R.Q, T.(Q), 16))
    end

    # the quotient against the resolvent, on every tabulated tableau
    tabs = Any[tab(BigFloat) for tab in TableauList]
    for s in 2:4,
        family in (TableauGauss, TableauRadauIA, TableauRadauIIA, TableauLobattoIIIA,
            TableauLobattoIIIC, TableauLobattoIIIF, TableauLobattoIIIG)

        push!(tabs, family(BigFloat, s))
    end
    for tab in tabs,
        z in (complex(BigFloat(-3) / 10, BigFloat(7) / 10), complex(BigFloat(-5) / 2))

        R = stability_function(tab)
        ref = resolvent(tab, z)
        @test abs(R(z) - ref) <= 1e-60 * abs(ref)
    end

    # symbolic coefficients
    R = stability_function(TableauGauss(symtype(), 2))
    @test R isa StabilityFunction{symtype()}
    @test all(iszero, SymPyPythonCall.simplify.(R.P .- [1, 1 // 2, 1 // 12]))
    @test all(iszero, SymPyPythonCall.simplify.(R.Q .- [1, -1 // 2, 1 // 12]))

    # a coefficient that is not finite
    for tab in (Tableau(:bad, 1, [NaN;;], [1.0], [0.0]), Tableau(
        :bad, 1, [0.0;;], [Inf], [0.0]))
        @test_throws r"^ArgumentError: the tableau has a coefficient that is not finite$" stability_function(tab)
        for f in (linear_order, isalgebraicallystable, imaginary_stability_interval ∘
                                                       PartitionedTableau)
            @test_throws r"^ArgumentError: the tableau has a coefficient that is not finite$" f(tab)
        end
        for pair in (PartitionedTableau(:x, tab, TableauExplicitEuler()),
            PartitionedTableau(:x, TableauExplicitEuler(), tab))
            @test_throws r"^ArgumentError: the tableau has a coefficient that is not finite$" imaginary_stability_interval(pair)
        end
    end
end

@testset "$(rpad("R(∞), the linear order and the Padé type",80))" begin
    # R(∞) of the closed forms
    for (tab, R∞) in ((TableauExplicitEuler(RB), 1 // 0), (TableauRK4(RB), 1 // 0),
        (TableauImplicitEuler(RB), 0), (TableauImplicitMidpoint(RB), -1),
        (TableauCrankNicolson(RB), -1), (zerotableau(RB), 1))
        r = stability_at_infinity(tab)
        @test r isa RB
        @test r == R∞
    end

    # R(∞) of the families, with round-off in the zero coefficients
    for T in (BigFloat, Float64)
        tol(x) = T == BigFloat ? 1e-70 : 16 * eps(Float64) * abs(x)
        for s in 1:4
            r = stability_at_infinity(TableauGauss(T, s))
            @test r isa T
            @test abs(r - (-1)^s) <= tol((-1)^s)
        end
        for s in 2:4
            for (family, R∞) in ((TableauRadauIB, (-1)^s), (
                TableauLobattoIIIA, (-1)^(s - 1)))
                r = stability_at_infinity(family(T, s))
                @test r isa T
                @test abs(r - R∞) <= tol(R∞)
            end
            for family in (TableauRadauIA, TableauRadauIIA, TableauLobattoIIIC)
                r = stability_at_infinity(family(T, s))
                @test r isa T
                @test iszero(r)
            end
            r = stability_at_infinity(TableauLobattoIII(T, s))
            @test r isa T
            @test r == T(Inf)
        end
    end

    # a small coefficient that is not round-off
    for θ in (1e-8, 1e8)
        r = stability_at_infinity(θmethod(Float64, θ))
        @test r isa Float64
        @test within(r, Float64(-(1 - big(θ)) / big(θ)), 16)
    end
    # in a rational type nothing counts as zero that is not 0
    θ = 1 // big(10)^30
    @test stability_at_infinity(θmethod(RB, θ)) == -(1 - θ) / θ

    # the linear order and the Padé type
    for T in (BigFloat, Float64), s in 1:4

        @test linear_order(TableauGauss(T, s)) === 2s
        @test pade_type(TableauGauss(T, s)) === (s, s)
    end
    for s in 2:4
        for (family, q, type) in ((TableauRadauIB, 2s, (s, s)), (
            TableauLobattoIIIF, 2s, (s, s)),
            (TableauRadauIA, 2s - 1, (s - 1, s)), (TableauRadauIIA, 2s - 1, (s - 1, s)),
            (TableauLobattoIIIA, 2s - 2, (s - 1, s - 1)),
            (TableauLobattoIIIC, 2s - 2, (s - 2, s)), (TableauLobattoIIIG, 2s - 2, nothing))
            tab = family(BigFloat, s)
            @test linear_order(tab) === q
            @test pade_type(tab) === type
        end
    end
    for (tab, q, type) in ((TableauExplicitEuler(RB), 1, (1, 0)),
        (TableauImplicitEuler(RB), 1, (0, 1)), (TableauKutta(RB), 3, (3, 0)),
        (TableauRK4(RB), 4, (4, 0)), (TableauRK5(RB), 5, nothing), (
        zerotableau(RB), 0, (0, 0)),
        (inconsistent(RB), 0, nothing))
        @test linear_order(tab) === q
        @test pade_type(tab) === type
    end
end

@testset "$(rpad("A- and L-stability",80))" begin
    # A-stable and L-stable tableaus
    for (T, s) in ((BigFloat, 1), (BigFloat, 2), (BigFloat, 3), (BigFloat, 4), (Float64, 1),
        (Float64, 2), (Float64, 3), (Float64, 4))
        @test isastable(TableauGauss(T, s))
        @test !islstable(TableauGauss(T, s))
    end
    for s in 2:4
        for family in (TableauRadauIA, TableauRadauIIA, TableauLobattoIIIC)
            @test isastable(family(BigFloat, s))
            @test islstable(family(BigFloat, s))
        end
        @test isastable(TableauLobattoIIIA(BigFloat, s))
        @test !islstable(TableauLobattoIIIA(BigFloat, s))
    end
    for tab in (TableauImplicitEuler(RB), θmethod(RB, 1))
        @test isastable(tab)
        @test islstable(tab)
    end
    for tab in (TableauImplicitMidpoint(RB), TableauCrankNicolson(RB), θmethod(RB, 1 // 2))
        @test isastable(tab)
        @test !islstable(tab)
    end
    @test isastable(zerotableau(RB))

    # tableaus that are not A-stable, one condition each
    for tab in (TableauExplicitEuler(RB), TableauRK4(RB), θmethod(RB, 1 // 4), poletableau(RB),
        TableauLobattoIII(BigFloat, 2), TableauLobattoIII(BigFloat, 3),
        TableauLobattoIII(BigFloat, 4))
        @test !isastable(tab)
        @test !islstable(tab)
    end
    # Q(z) = P(z) = 1 + z², E ≡ 0: only the Routh test sees the poles ±i on the axis, where its
    # first column holds a zero
    @test !isastable(Tableau{RB}(:axis, 1, [0 1; -1 0], [0, 0], [0, 0]))

    # a symbolic tableau
    tab = TableauGauss(symtype(), 2)
    for f in (isastable, islstable, linear_order, real_stability_interval)
        @test_throws Regex("^MethodError: no method matching $(nameof(f))\\(") f(tab)
    end
end

@testset "$(rpad("The stability intervals",80))" begin
    # closed forms on the real axis
    kutta = newton_root([-12, 6, -3, 1], 2.5)
    rk4 = newton_root([-24, 12, -4, 1], 2.8)
    for T in (RB, BigFloat, Float64)
        F = float(T)
        for (tab, r) in ((TableauExplicitEuler(T), 2), (TableauKutta(T), kutta), (
            TableauRK4(T), rk4))
            x = real_stability_interval(tab)
            @test x isa F
            @test within(x, F(r), 16)
        end
    end
    for T in (BigFloat, Float64)
        x = real_stability_interval(TableauLobattoIII(T, 2))
        @test x isa T
        @test within(x, T(2), 16)
    end
    for s in 1:4
        @test isvalue(real_stability_interval(TableauGauss(BigFloat, s)), BigFloat(Inf))
    end
    for s in 2:4
        @test isvalue(real_stability_interval(TableauRadauIIA(BigFloat, s)), BigFloat(Inf))
    end

    # a touching point inside the real interval
    for T in (RB, Float64)
        x = real_stability_interval(rkc2(T))
        @test x isa float(T)
        @test within(x, float(T)(8), 16)
    end

    # closed forms on the imaginary axis
    for T in (RB, Float64)
        F = float(T)
        for tab in (TableauExplicitEuler(T), TableauRK5(T))
            y = imaginary_stability_interval(tab)
            @test y isa F
            @test isvalue(y, zero(F))
        end
        for (tab, r) in ((TableauKutta(T), sqrt(big(3))), (TableauRK4(T), 2 * sqrt(big(2))))
            y = imaginary_stability_interval(tab)
            @test y isa F
            @test within(y, F(r), 16)
        end
    end
    @test isvalue(imaginary_stability_interval(TableauLobattoIII(BigFloat, 3)), zero(BigFloat))
    for s in 2:4, family in (TableauGauss, TableauRadauIIA, TableauLobattoIIIC)

        @test isvalue(imaginary_stability_interval(family(BigFloat, s)), BigFloat(Inf))
    end

    # the zero tableau: R ≡ 1
    @test isvalue(real_stability_interval(zerotableau(RB)), BigFloat(Inf))
    @test isvalue(imaginary_stability_interval(zerotableau(RB)), BigFloat(Inf))

    # the precision of the result
    @test real_stability_interval(TableauGauss(Float32, 2)) === Inf32
    @test imaginary_stability_interval(TableauGauss(Float32, 2)) === Inf32
    @test real_stability_interval(TableauExplicitEuler(Float32)) === 2.0f0
    @test imaginary_stability_interval(TableauExplicitEuler(Float32)) === 0.0f0
end

@testset "$(rpad("The stability interval of a partitioned tableau",80))" begin
    # a tableau paired with itself
    for (tab, r) in ((TableauExplicitEuler(RB), big(0)), (TableauRK4(RB), 2 * sqrt(big(2))),
        (TableauGauss(BigFloat, 2), BigFloat(Inf)))
        y = imaginary_stability_interval(PartitionedTableau(tab))
        y1 = imaginary_stability_interval(tab)
        @test y isa BigFloat
        if isfinite(r) && !iszero(r)
            @test within(y, r, 16)
            @test within(y, y1, 16)
        else
            @test y == r
            @test y == y1
        end
    end

    # the symplectic Euler method and the Störmer–Verlet pair
    for T in (RB, Float64)
        tab = PartitionedTableau(:SymplecticEuler, TableauExplicitEuler(T), TableauImplicitEuler(T))
        y = imaginary_stability_interval(tab)
        @test y isa float(T)
        @test within(y, float(T)(2), 16)
    end
    y = imaginary_stability_interval(TableauLobattoIIIAIIIB(BigFloat, 2))
    @test y isa BigFloat
    @test within(y, BigFloat(2), 16)

    # M(ν) = [1 ν; ν 1], with det M = 1 - ν² and the eigenvalue 1 + ν > 1 for every ν > 0: only
    # the condition det(I - M) ≥ 0 ends the interval, at 0
    tab = PartitionedTableau(:hyperbolic, TableauExplicitEuler(RB),
        Tableau{RB}(:minus, 1, [0;;], [-1], [0]))
    @test isvalue(imaginary_stability_interval(tab), zero(BigFloat))

    # two roots closer than the precision of Float64 are both found
    roots = RungeKutta._positive_roots(RB[1 + big(2) ^ -60, -2 - big(2) ^ -60, 1], Float64)
    @test length(roots) == 2
    @test all(within.(Float64.(roots), 1.0, 1))

    # every other name raises a MethodError for a pair
    pair = PartitionedTableau(TableauGauss(2))
    for f in (stability_function, stability_at_infinity, linear_order, pade_type, isastable,
        islstable, isalgebraicallystable, algebraic_stability_matrix, real_stability_interval,
        stability_region)
        @test_throws Regex("^MethodError: no method matching $(nameof(f))\\(") f(pair)
    end
end

@testset "$(rpad("Algebraic stability",80))" begin
    # the matrix
    for θ in (1 // 4, 1 // 2, 1 // 1)
        M = algebraic_stability_matrix(θmethod(RB, θ))
        @test M isa Matrix{RB}
        @test M == [-(1 - θ)^2 0; 0 θ^2]
    end
    for s in 1:4
        M = algebraic_stability_matrix(TableauGauss(BigFloat, s))
        @test M isa Matrix{BigFloat}
        @test all(abs.(M) .<= 1e-70)
    end

    # the predicate
    for s in 1:4
        @test isalgebraicallystable(TableauGauss(BigFloat, s))
    end
    for s in 2:4
        for family in (TableauRadauIA, TableauRadauIIA, TableauLobattoIIIC)
            @test isalgebraicallystable(family(BigFloat, s))
        end
        @test !isalgebraicallystable(TableauLobattoIIIA(BigFloat, s))
    end
    for tab in (TableauImplicitEuler(RB), θmethod(RB, 1))
        @test isalgebraicallystable(tab)
    end
    for tab in (θmethod(RB, 1 // 2), TableauExplicitEuler(RB), TableauRK4(RB), poletableau(RB))
        @test !isalgebraicallystable(tab)
    end
end

@testset "$(rpad("The data of the stability region",80))" begin
    # a grid through a pole
    for T in (Float64, RB)
        F = float(T)
        x, y, absR = stability_region(TableauImplicitEuler(T); xlims = (-1.0, 3.0),
            ylims = (-2.0, 2.0), n = 5)
        @test x isa Vector{F}
        @test y isa Vector{F}
        @test absR isa Matrix{F}
        @test x == [-1, 0, 1, 2, 3]
        @test y == [-2, -1, 0, 1, 2]
        @test absR[3, 3] == F(Inf)
        @test absR[2, 3] == 1
        @test absR[1, 3] == 1 // 2
    end

    # the default window and the values
    tab = TableauRK4()
    x, y, absR = stability_region(tab)
    @test length(x) == length(y) == 201
    @test x[1] == -6
    @test x[end] == 2
    @test y[1] == -4
    @test y[end] == 4
    @test size(absR) == (201, 201)
    for i in 1:50:201, j in 1:50:201

        ref = Float64(abs(resolvent(TableauRK4(BigFloat), complex(big(x[i]), big(y[j])))))
        @test within(absR[i, j], ref, 16)
    end

    # a window that is refused
    for kw in ((; xlims = (1.0, 1.0)), (; xlims = (2.0, 1.0)), (; ylims = (NaN, 1.0)),
        (; ylims = (-Inf, 0.0)), (; n = 1))
        @test_throws r"^ArgumentError: the window must have finite, increasing limits and n ≥ 2$" stability_region(
            tab; kw...)
    end
end
