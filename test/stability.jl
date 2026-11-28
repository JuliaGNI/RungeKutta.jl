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

    # R(z) at a real and a complex z in Float64, and at a rational z, each in the type of z
    R = stability_function(TableauRK4(Float64))
    for z in (-0.7, 2.3, complex(-0.3, 0.7), complex(-2.5, 1.5))
        ref = resolvent(TableauRK4(BigFloat), big(z))
        @test R(z) isa typeof(z)
        @test within(real(R(z)), Float64(real(ref)), 16)
        @test within(imag(R(z)), Float64(imag(ref)), 16)
    end
    R = stability_function(TableauImplicitMidpoint(RB))
    for z in (3 // big(5), -7 // big(2))
        @test R(z) isa RB
        @test R(z) == (1 + z / 2) / (1 - z / 2)
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
    # R(z) = 1 + a z + a² z², E(y) = a² y² (1 - a² y²): E < 0 for y > 1/a = 2^2048, a root
    # beyond floatmax(Float64), which an exact decision does not convert
    a = 1 // big(2)^2048
    @test isastable(Tableau{RB}(:far, 2, [0 0; a 0], [0, a], [0, 0])) === false
    # R(z) = (1 - z/2 - z²/4) / ((1 - z)(1 - z/2)), E(y) = y² (8 + 3y²) / 16 has no positive
    # root: the Sturm counts of 8 + 3y² at 0 and at the Cauchy bound agree
    diag = Tableau{RB}(:diag, 1, [1 0; 0 1 // 2], [1 // 2, 1 // 2], [1, 1 // 2])
    @test isastable(diag) === true
    @test islstable(diag) === false

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
    # R(x) = T₃(1 + x/9) touches |R| = 1 at x = -4.5 and -13.5, where the rounding of b in a
    # floating-point type splits each touching point into two close roots; the endpoint -18 has
    # the condition number 272
    cheb(T) = Tableau{T}(
        :cheb, 1, [0 0 0; 1 0 0; 0 1 0], [621 // 729, 104 // 729, 4 // 729],
        [0, 1, 1])
    for T in (RB, Float64, Float32)
        x = real_stability_interval(cheb(T))
        @test x isa float(T)
        @test within(x, float(T)(18), 1024)
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
    # an explicit tableau with R(z) = 1 + z + βz² + γz³ + z⁴/10, whose
    # E(y) = -y²(y² - 1)²(y² - u₀)/100 touches zero at y = 1; the rounding of its coefficients
    # splits the touching point into two close roots, at which the interval must not end
    β = big"0.5262500227059282687825755391373232433684"
    γ = big"0.1809694977871384716612685727462138383356"
    u₀ = big"5.2500045411856537565151078274646486736725"
    for T in (Float64, Float32)
        b, c, d = T(β), T(γ), T(1 // 10)
        tab = Tableau{T}(:touch, 1, [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0],
            [1 - b, b - c, c - d, d], [0, 1, 1, 1])
        exact = Tableau{RB}(
            :touch, 1, RB.(Matrix(tab.a)), RB.(Vector(tab.b)), RB.(Vector(tab.c)))
        @test imaginary_stability_interval(exact) < 1.01
        y = imaginary_stability_interval(tab)
        @test y isa T
        @test within(y, T(sqrt(u₀)), 16)
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

    # Gauss(4) paired with itself in Float64: no round-off of D, N1 or N2 ends the interval
    @test imaginary_stability_interval(PartitionedTableau(TableauGauss(Float64, 4))) === Inf

    # every other name raises a MethodError for a pair
    pair = PartitionedTableau(TableauGauss(2))
    for f in (stability_function, stability_at_infinity, linear_order, pade_type, isastable,
        islstable, isalgebraicallystable, algebraic_stability_matrix, real_stability_interval,
        stability_region)
        @test_throws Regex("^MethodError: no method matching $(nameof(f))\\(") f(pair)
    end
end

@testset "$(rpad("Float32 decisions",80))" begin
    F = Float32
    # A-stability
    for s in 1:4
        @test isastable(TableauGauss(F, s))
    end
    for s in 2:4
        for family in (TableauRadauIA, TableauRadauIIA, TableauLobattoIIIA, TableauLobattoIIIC)
            @test isastable(family(F, s))
        end
        @test !isastable(TableauLobattoIII(F, s))
    end
    # the constant coefficient Q(0) = 1 never counts as zero, at a scale beyond 1/rtol
    @test isastable(θmethod(Float64, 1e15))

    # the intervals of the single tableaus
    gauss = [TableauGauss(F, s) for s in 1:4]
    radau = [TableauRadauIIA(F, s) for s in 2:4]
    lobatto = [TableauLobattoIIIC(F, s) for s in 2:4]
    for tab in vcat(gauss, radau, lobatto)
        @test imaginary_stability_interval(tab) === Inf32
    end
    for tab in vcat(gauss, radau)
        @test real_stability_interval(tab) === Inf32
    end
    for (tab, r) in ((TableauRK4(F), 2 * sqrt(big(2))), (TableauKutta(F), sqrt(big(3))))
        y = imaginary_stability_interval(tab)
        @test y isa F
        @test within(y, F(r), 16)
    end
    for (tab, c, x0) in ((TableauRK4(F), [-24, 12, -4, 1], 2.8),
        (TableauKutta(F), [-12, 6, -3, 1], 2.5))
        x = real_stability_interval(tab)
        @test x isa F
        @test within(x, F(newton_root(c, x0)), 16)
    end
    # the endpoint of TableauRK5 has the condition number 778
    x = real_stability_interval(TableauRK5(F))
    @test x isa F
    @test within(x, F(real_stability_interval(TableauRK5(RB))), 1024)

    # a tableau paired with itself
    for tab in vcat(
        [TableauExplicitEuler(F), TableauKutta(F), TableauRK4(F)], gauss, radau,
        lobatto)
        y = imaginary_stability_interval(PartitionedTableau(tab))
        y1 = imaginary_stability_interval(tab)
        @test y isa F
        if isfinite(y1) && !iszero(y1)
            @test within(y, y1, 16)
        else
            @test y === y1
        end
    end
end

@testset "$(rpad("The round-off bounds of the decisions",80))" begin
    B = RungeKutta._BoundedPoly
    r = 1 // 100

    # P and Q of the θ-method, θ = 1/4: the constant and a zero coefficient carry no bound, every
    # other coefficient rtol · m, with m = 1
    rtol = 16 * eps(Float64)
    P, Q = RungeKutta._bounded_stability_function(θmethod(Float64, 0.25), rtol)
    @test P.v == [1, 0.75, 0] && Q.v == [1, -0.25, 0]
    @test P.e == [0, rtol, 0] && Q.e == [0, rtol, 0]

    # a product: Σ (aₓᵢ e_yⱼ + eₓᵢ a_yⱼ + rtol aₓᵢ a_yⱼ) over i + j = k, in exact arithmetic
    x = B(RB[1, 2], RB[1, 2], RB[0, 1 // 10])
    y = B(RB[3, -1], RB[3, 1], RB[1 // 10, 1 // 5])
    xy = RungeKutta._bmul(x, y, r)
    @test xy.v == [3, 5, -2] && xy.a == [3, 7, 2]
    @test xy.e == [13 // 100, 77 // 100, 13 // 25]
    # a negative factor scales the bounds by its magnitude
    @test RungeKutta._bscale(-2, x).e == [0, 1 // 5]

    # Faddeev–LeVerrier with its bounds on diag(1, 2): det(I - z M) = 1 - 3z + 2z², with the
    # bounds 3 rtol (|a| + |d|) and (5 rtol (a² + d²) + 14 rtol |a d|) / 2
    D = RungeKutta._det_coefficients_bounded(RB[1 0; 0 2], r)
    @test D.v == [1, -3, 2]
    @test D.e == [0, 9 * r, 53 * r / 2]

    # |p(iy)|² is even in y: its odd coefficients are exactly zero, also where the product
    # p(z) p(-z) holds round-off in them
    for tab in TableauList, p in RungeKutta._bounded_stability_function(tab(Float64), rtol)

        E = RungeKutta._abs2_on_imaginary_axis(p, rtol)
        @test all(iszero, E.v[2:2:end]) && all(iszero, E.e[2:2:end])
    end

    # roots closer than the bisection can separate are an error, not one root
    ε = 1 // big(2)^5000
    @test_throws r"^the bisection did not isolate a root of the Sturm sequence$" RungeKutta._positive_roots(
        RB[1 + ε, -2 - ε, 1], Float64)
    # an endpoint 2^-10000 that the bisection cannot reach in its step budget is an error, not
    # a wrong value
    @test_throws r"^the bisection did not locate a root to the precision of BigFloat" real_stability_interval(Tableau{
        RB}(
        :k19, 1, RB[0;;], RB[big(2) ^ 10001], RB[0]))
end

@testset "$(rpad("Algebraic stability",80))" begin
    # the matrix
    for θ in (1 // 4, 1 // 2, 1 // 1)
        M = algebraic_stability_matrix(θmethod(RB, θ))
        @test M isa Matrix{RB}
        @test M == [-(1 - θ)^2 0; 0 θ^2]
    end
    for θ in (0.25, 0.5, 1.0)
        M = algebraic_stability_matrix(θmethod(Float64, θ))
        @test M isa Matrix{Float64}
        @test all(within.(M, [-(1 - θ)^2 0; 0 θ^2], 16))
    end
    for T in (BigFloat, Float64), s in 1:4

        M = algebraic_stability_matrix(TableauGauss(T, s))
        @test M isa Matrix{T}
        @test all(abs.(M) .<= (T == BigFloat ? 1e-70 : 16 * eps(Float64)))
    end

    # the predicate
    for T in (BigFloat, Float64)
        for s in 1:4
            @test isalgebraicallystable(TableauGauss(T, s))
        end
        for s in 2:4
            for family in (TableauRadauIA, TableauRadauIIA, TableauLobattoIIIC)
                @test isalgebraicallystable(family(T, s))
            end
            @test !isalgebraicallystable(TableauLobattoIIIA(T, s))
        end
    end
    for T in (RB, Float64)
        for tab in (TableauImplicitEuler(T), θmethod(T, 1))
            @test isalgebraicallystable(tab)
        end
        for tab in (θmethod(T, 1 // 2), TableauExplicitEuler(T), TableauRK4(T), poletableau(T))
            @test !isalgebraicallystable(tab)
        end
    end
    # M = (1f20 - 1/4) [1 1; 1 1], positive semidefinite of rank one: the elimination divides
    # the pivot column by the pivot before the product, whose entries 1f40 overflow Float32
    k23 = Tableau{Float32}(:k23, 1, fill(1.0f20, 2, 2), Float32[1 / 2, 1 / 2],
        Float32[2.0f20, 2.0f20])
    @test isalgebraicallystable(k23) === true

    # the inferred return type
    for T in (Float64, Float32, BigFloat, RB)
        @test @inferred(isalgebraicallystable(TableauGauss(T, 2))) isa Bool
        @test @inferred(isalgebraicallystable(TableauRK4(T))) isa Bool
    end
end

@testset "$(rpad("The data of the stability region",80))" begin
    # a grid through a pole
    for T in (Float64, BigFloat, RB)
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
    # the limits are checked in the precision of the grid, float(T)
    for kw in ((; xlims = (-1e300, 1e300)), (; ylims = (1.0, 1.0 + 1e-12)))
        @test_throws r"^ArgumentError: the window must have finite, increasing limits and n ≥ 2$" stability_region(
            TableauRK4(Float32); kw...)
    end
end

# A variable that a closure captures and the enclosing method reassigns is a `Core.Box`, and every
# operation on it infers `Any`. The lowered code of each method of `src/stability.jl` holds none.
@testset "$(rpad("No captured variable is boxed",80))" begin
    scanned = Symbol[]
    boxed = Tuple{Symbol, Int}[]
    for name in names(RungeKutta; all = true)
        isdefined(RungeKutta, name) || continue
        f = getfield(RungeKutta, name)
        f isa Function || continue
        for m in methods(f)
            (m.module === RungeKutta && endswith(string(m.file), "stability.jl")) ||
                continue
            push!(scanned, name)
            occursin("Core.Box", string(Base.uncompressed_ir(m))) &&
                push!(boxed, (m.name, Int(m.line)))
        end
    end
    @test isempty(boxed)
    # the scan reaches the methods of the file
    @test :_has_roots_in_right_half_plane_only in scanned
    @test :isalgebraicallystable in scanned
    @test length(scanned) >= 20
end
