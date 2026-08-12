using QuadratureRules: radau_legendre_nodes, radau_legendre_weights
using RungeKutta.Tableaus: radau_1_coefficients, radau_2_coefficients

@testset "$(rpad("Radau Tableaus",80))" begin

    # The one-node Radau quadrature rule is perfectly well defined — it is a Riemann sum —
    # so the nodes and weights do not throw. What is undefined is the one-stage Radau
    # *tableau*, and that restriction lives in the coefficients and hence in the constructors.
    @test_throws ErrorException radau_1_coefficients(1)
    @test_throws ErrorException radau_2_coefficients(1)

    @test_throws ErrorException TableauRadauIA(1)
    @test_throws ErrorException TableauRadauIB(1)
    @test_throws ErrorException TableauRadauIIA(1)
    @test_throws ErrorException TableauRadauIIB(1)

    
    function _TableauRadauIA2(T=Float64)
        a = [[1//4  -1//4  ]
             [1//4   5//12 ]]
        b = [1//4, 3//4]
        c = [0,    2//3]

        Tableau{T}(:RadauIA2, 3, a, b, c; R∞=0)
    end

    function _TableauRadauIA3(T=Float64)
        a = [
                [ 1//9     (- 1 -    √6)/18    (- 1 +    √6)/18  ]
                [ 1//9     ( 88 +  7*√6)/360   ( 88 - 43*√6)/360 ]
                [ 1//9     ( 88 + 43*√6)/360   ( 88 -  7*√6)/360 ]
            ]
        b = [1/9,         (16+√6)/36,  (16-√6)/36 ]
        c = [0,           ( 6-√6)/10,  ( 6+√6)/10 ]

        Tableau{T}(:RadauIA3, 5, a, b, c; R∞=0)
    end

    function _TableauRadauIB2(T=Float64)
        a = [[1//8  -1//8  ]
             [7//24  3//8  ]]
        b = [1//4, 3//4]
        c = [0,    2//3]

        Tableau{T}(:RadauIB2, 3, a, b, c; R∞=0)
    end

    function _TableauRadauIB3(T=Float64)
        a = [
                [            1//18    (- 1 -     √6)/36    (- 1 +     √6)/36   ]
                [ ( 52 +  3*√6)/450   ( 16 +     √6)/72    (472 - 217*√6)/1800 ]
                [ ( 52 -  3*√6)/450   (472 + 217*√6)/1800  ( 16 -     √6)/72   ]
            ]
        b = [1/9,         (16+√6)/36,  (16-√6)/36 ]
        c = [0,           ( 6-√6)/10,  ( 6+√6)/10 ]

        Tableau{T}(:RadauIB3, 5, a, b, c; R∞=0)
    end

    function _TableauRadauIIA2(T=Float64)
        a = [[5//12  -1//12]
             [3//4    1//4 ]]
        b = [3//4, 1//4]
        c = [1//3, 1//1]

        Tableau{T}(:RadauIIA2, 3, a, b, c; R∞=0)
    end

    function _TableauRadauIIA3(T=Float64)
        a = [
                [ 11/45 -  7*√6/360    37/225-169*√6/1800   -2/225+√6/75 ]
                [ 37/225+169*√6/1800   11/45 +  7*√6/360    -2/225-√6/75 ]
                [  4/9  -    √6/36      4/9  +    √6/36      1/9         ]
            ]
        b = [4/9-√6/36,   4/9+√6/36,  1/9 ]
        c = [2/5-√6/10,   2/5+√6/10,  1   ]

        Tableau{T}(:RadauIIA3, 5, a, b, c; R∞=0)
    end

    function _TableauRadauIIB2(T=Float64)
        a = [[3//8   -1//24]
             [7//8    1//8 ]]
        b = [3//4, 1//4]
        c = [1//3, 1//1]

        Tableau{T}(:RadauIIB2, 3, a, b, c; R∞=0)
    end

    function _TableauRadauIIB3(T=Float64)
        a = [
                [( 16 -     √6)/72    (328 - 167*√6)/1800   (-2 + 3*√6)/450 ]
                [(328 + 167*√6)/1800  ( 16 +     √6)/72     (-2 - 3*√6)/450 ]
                [( 85 -  10*√6)/180   ( 85 +  10*√6)/180              1/18  ]
            ]
        b = [4/9-√6/36,   4/9+√6/36,  1/9 ]
        c = [2/5-√6/10,   2/5+√6/10,  1   ]

        Tableau{T}(:RadauIIB3, 5, a, b, c; R∞=0)
    end

    @test_throws ErrorException TableauRadauIA(1)
    @test_throws ErrorException TableauRadauIIA(1)

    @test TableauRadauIA(2)  ≈ _TableauRadauIA2()
    @test TableauRadauIA(3)  ≈ _TableauRadauIA3()

    @test TableauRadauIB(2)  ≈ _TableauRadauIB2()
    @test TableauRadauIB(3)  ≈ _TableauRadauIB3()

    @test TableauRadauIIA(2) ≈ _TableauRadauIIA2()
    @test TableauRadauIIA(3) ≈ _TableauRadauIIA3()

    @test TableauRadauIIB(2) ≈ _TableauRadauIIB2()
    @test TableauRadauIIB(3) ≈ _TableauRadauIIB3()

    @test reference(TableauRadauIA(2))  == reference(Val(:RadauIA))
    @test reference(TableauRadauIB(2))  == reference(Val(:RadauIB))
    @test reference(TableauRadauIIA(2)) == reference(Val(:RadauIIA))
    @test reference(TableauRadauIIB(2)) == reference(Val(:RadauIIB))

    for s in 2:5
        @test !issymplectic(TableauRadauIA(s))
        @test !issymplectic(TableauRadauIIA(s))
    end

    for T in (Float32, Float64, BigFloat, symtype())
        for s in 2:3
            @test_nowarn radau_legendre_nodes(T, s, Val(:left))
            @test_nowarn radau_legendre_weights(T, s, Val(:left))
            @test_nowarn radau_1_coefficients(T,s)

            @test_nowarn radau_legendre_nodes(T, s, Val(:right))
            @test_nowarn radau_legendre_weights(T, s, Val(:right))
            @test_nowarn radau_2_coefficients(T,s)

            @test_nowarn TableauRadauIA(T,s)
            @test_nowarn TableauRadauIIA(T,s)
        end
    end

    @test radau_legendre_nodes(Float32, 2, Val(:left)) ≈ radau_legendre_nodes(Float64, 2, Val(:left))
    @test radau_legendre_weights(Float32, 2, Val(:left)) ≈ radau_legendre_weights(Float64, 2, Val(:left))
    @test radau_1_coefficients(Float32,2) ≈ radau_1_coefficients(Float64,2)

    @test radau_legendre_nodes(Float32, 2, Val(:right)) ≈ radau_legendre_nodes(Float64, 2, Val(:right))
    @test radau_legendre_weights(Float32, 2, Val(:right)) ≈ radau_legendre_weights(Float64, 2, Val(:right))
    @test radau_2_coefficients(Float32,2) ≈ radau_2_coefficients(Float64,2)

    @test radau_legendre_nodes(symtype(), 2, Val(:left)) ≈ radau_legendre_nodes(Float64, 2, Val(:left))
    @test radau_legendre_weights(symtype(), 2, Val(:left)) ≈ radau_legendre_weights(Float64, 2, Val(:left))
    @test radau_1_coefficients(symtype(),2) ≈ radau_1_coefficients(Float64,2)

    @test radau_legendre_nodes(symtype(), 2, Val(:right)) ≈ radau_legendre_nodes(Float64, 2, Val(:right))
    @test radau_legendre_weights(symtype(), 2, Val(:right)) ≈ radau_legendre_weights(Float64, 2, Val(:right))
    @test radau_2_coefficients(symtype(),2) ≈ radau_2_coefficients(Float64,2)

    @test TableauRadauIA(Float32,2) ≈ TableauRadauIA(Float64,2)
    @test TableauRadauIIA(Float32,2) ≈ TableauRadauIIA(Float64,2)
    @test TableauRadauIA(symtype(),2) ≈ TableauRadauIA(Float64,2)
    @test TableauRadauIIA(symtype(),2) ≈ TableauRadauIIA(Float64,2)

    # The s-stage Radau quadrature integrates polynomials up to degree 2s-2
    # exactly. This checks the arbitrary precision nodes and weights directly,
    # rather than only via their double precision counterparts. The Radau IA
    # nodes include the left endpoint, the Radau IIA nodes the right one.
    for s in 2:10
        b₁ = radau_legendre_weights(BigFloat, s, Val(:left))
        c₁ = radau_legendre_nodes(BigFloat, s, Val(:left))
        b₂ = radau_legendre_weights(BigFloat, s, Val(:right))
        c₂ = radau_legendre_nodes(BigFloat, s, Val(:right))

        @test eltype(b₁) == eltype(c₁) == eltype(b₂) == eltype(c₂) == BigFloat
        @test issorted(c₁) && issorted(c₂)
        @test c₁[begin] == 0 && c₁[end] < 1
        @test c₂[begin] > 0 && c₂[end] == 1
        @test c₁ ≈ 1 .- reverse(c₂)

        for k in 0:2s-2
            @test sum(b₁ .* c₁.^k) ≈ 1 / BigFloat(k+1) atol=1E-60
            @test sum(b₂ .* c₂.^k) ≈ 1 / BigFloat(k+1) atol=1E-60
        end
    end

end
