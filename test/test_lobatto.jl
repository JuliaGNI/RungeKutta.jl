import LinearAlgebra

using RungeKutta.Tableaus: get_lobatto_nodes, get_lobatto_weights,
                  get_lobatto_a_coefficients, get_lobatto_b_coefficients,
                  get_lobatto_c_coefficients, get_lobatto_c̄_coefficients,
                  get_lobatto_d_coefficients, get_lobatto_e_coefficients,
                  get_lobatto_f_coefficients, get_lobatto_g_coefficients

@testset "$(rpad("Lobatto Tableaus",80))" begin

    @test_throws ErrorException get_lobatto_nodes(1)
    @test_throws ErrorException get_lobatto_weights(1)
    @test_throws ErrorException get_lobatto_nullvector(1)
    @test_throws ErrorException get_lobatto_a_coefficients(1)
    @test_throws ErrorException get_lobatto_b_coefficients(1)
    @test_throws ErrorException get_lobatto_c_coefficients(1)
    @test_throws ErrorException get_lobatto_c̄_coefficients(1)
    @test_throws ErrorException get_lobatto_d_coefficients(1)
    @test_throws ErrorException get_lobatto_e_coefficients(1)
    @test_throws ErrorException get_lobatto_f_coefficients(1)
    @test_throws ErrorException get_lobatto_g_coefficients(1)


    function _get_lobatto_nullvector(s)
        if s == 2
            d = [+1.0, -1.0]
        elseif s == 3
            d = [+1.0, -2.0, +1.0]
        elseif s == 4
            d = [+1.0, -√5, +√5, -1.0]
        elseif s == 5
            d = [+3.0, -7.0, +8.0, -7.0, +3.0]
        else
            @error("We don't have a d vector for s=$(s) stages.")
        end
        return LinearAlgebra.normalize(d) * sign(d[begin])
    end


    @test get_lobatto_nullvector(2; normalize=true) ≈ _get_lobatto_nullvector(2)
    @test get_lobatto_nullvector(3; normalize=true) ≈ _get_lobatto_nullvector(3)
    @test get_lobatto_nullvector(4; normalize=true) ≈ _get_lobatto_nullvector(4)
    @test get_lobatto_nullvector(5; normalize=true) ≈ _get_lobatto_nullvector(5)


    function _getTableauLobattoIIIA2(T=Float64)
        a = BigFloat[
                [0     0   ]
                [1/2   1/2 ]
            ]

        Tableau{T}(:LobattoIIIA2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2); R∞ = -1)
    end

    function _getTableauLobattoIIIA3(T=Float64)
        a = BigFloat[
                [0      0     0    ]
                [5/24   1/3  -1/24 ]
                [1/6    2/3   1/6  ]
            ]

        Tableau{T}(:LobattoIIIA3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3); R∞ = +1)
    end

    function _getTableauLobattoIIIA4(T=Float64)
        a = BigFloat[
                [0            0               0               0          ]
                [(11+√5)/120  (25-   √5)/120  (25-13*√5)/120  (-1+√5)/120]
                [(11-√5)/120  (25+13*√5)/120  (25+   √5)/120  (-1-√5)/120]
                [      1/12            5/12            5/12         1/12 ]
            ]

        Tableau{T}(:LobattoIIIA4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4); R∞ = -1)
    end


    function _getTableauLobattoIIIA5(T=Float64)
        a = BigFloat[
                [0                 0                   0                  0                   0               ]
                [(119+3*√21)/1960  (343-  9*√21)/2520  (392-96*√21)/2205  (343- 69*√21)/2520  (-21+3*√21)/1960]
                [         13/320   (392+105*√21)/2880             8/45    (392-105*√21)/2880            3/320 ]
                [(119-3*√21)/1960  (343+ 69*√21)/2520  (392+96*√21)/2205  (343+  9*√21)/2520  (-21-3*√21)/1960]
                [          1/20               49/180             16/45               49/180             1/20  ]
            ]

        Tableau{T}(:LobattoIIIA5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5); R∞ = +1)
    end


    function _getTableauLobattoIIIB2(T=Float64)
        a = BigFloat[
                [1/2  0]
                [1/2  0]
            ]

        Tableau{T}(:LobattoIIIB2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2); R∞ = -1)
    end

    function _getTableauLobattoIIIB3(T=Float64)
        a = BigFloat[
                [1/6  -1/6   0   ]
                [1/6   1/3   0   ]
                [1/6   5/6   0   ]
            ]

        Tableau{T}(:LobattoIIIB3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3); R∞ = +1)
    end

    function _getTableauLobattoIIIB4(T=Float64)
        a = BigFloat[
                [ 1/12  (-1-   √5)/24   (-1+   √5)/24    0 ]
                [ 1/12  (25+   √5)/120  (25-13*√5)/120   0 ]
                [ 1/12  (25+13*√5)/120  (25-   √5)/120   0 ]
                [ 1/12  (11-   √5)/24   (11+   √5)/24    0 ]
            ]

        Tableau{T}(:LobattoIIIB4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4); R∞ = -1)
    end


    function _getTableauLobattoIIIB5(T=Float64)
        a = BigFloat[
                [ 1/20  ( -7-   √21)/120             1/15   ( -7+   √21)/120    0 ]
                [ 1/20  (343+ 9*√21)/2520  (56-15*√21)/315  (343-69*√21)/2520   0 ]
                [ 1/20  ( 49+12*√21)/360             8/45   ( 49-12*√21)/360    0 ]
                [ 1/20  (343+69*√21)/2520  (56+15*√21)/315  (343- 9*√21)/2520   0 ]
                [ 1/20  (119- 3*√21)/360            13/45   (119+ 3*√21)/360    0 ]
            ]

        Tableau{T}(:LobattoIIIB5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5); R∞ = +1)
    end


    function _getTableauLobattoIIIC2(T=Float64)
        a = BigFloat[
                [1/2  -1/2 ]
                [1/2   1/2 ]
            ]

        Tableau{T}(:LobattoIIIC2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2); R∞ = -1)
    end

    function _getTableauLobattoIIIC3(T=Float64)
        a = BigFloat[
                [1/6  -1/3    1/6  ]
                [1/6   5/12  -1/12 ]
                [1/6   2/3    1/6  ]
            ]

        Tableau{T}(:LobattoIIIC3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3); R∞ = +1)
    end

    function _getTableauLobattoIIIC4(T=Float64)
        a = BigFloat[
                [ 1/12        -√5/12         √5/12   -1/12 ]
                [ 1/12          1/4   (10-7*√5)/60   √5/60 ]
                [ 1/12  (10+7*√5)/60          1/4   -√5/60 ]
                [ 1/12          5/12          5/12    1/12 ]
            ]

        Tableau{T}(:LobattoIIIC4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4); R∞ = -1)
    end


    function _getTableauLobattoIIIC5(T=Float64)
        a = BigFloat[
                [ 1/20            -21/180             2/15             -21/180     1/20  ]
                [ 1/20             29/180   (47-15*√21)/315  (203- 30*√21)/1260   -3/140 ]
                [ 1/20  (329+105*√21)/2880           73/360  (329-105*√21)/2880    3/160 ]
                [ 1/20  (203+ 30*√21)/1260  (47+15*√21)/315             29/180    -3/140 ]
                [ 1/20             49/180            16/45              49/180     1/20  ]
            ]

        Tableau{T}(:LobattoIIIC5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5); R∞ = +1)
    end


    function _getTableauLobattoIIIC̄2(T=Float64)
        a = BigFloat[
                [0  0]
                [1  0]
            ]

        Tableau{T}(:LobattoIIIC̄2, 2, a, get_lobatto_weights(2), get_lobatto_nodes(2); R∞ = -1)
    end

    function _getTableauLobattoIIIC̄3(T=Float64)
        a = BigFloat[
                [0    0    0 ]
                [1/4  1/4  0 ]
                [0    1    0 ]
            ]

        Tableau{T}(:LobattoIIIC̄3, 4, a, get_lobatto_weights(3), get_lobatto_nodes(3); R∞ = +1)
    end

    function _getTableauLobattoIIIC̄4(T=Float64)
        a = BigFloat[
                [      0             0             0     0 ]
                [ (5+√5)/60          1/6   (15-7*√5)/60  0 ]
                [ (5-√5)/60  (15+7*√5)/60          1/6   0 ]
                [      1/6      (5-√5)/12     (5+√5)/12  0 ]
            ]

        Tableau{T}(:LobattoIIIC̄4, 6, a, get_lobatto_weights(4), get_lobatto_nodes(4); R∞ = -1)
    end


    function _getTableauLobattoIIIC̄5(T=Float64)
        a = BigFloat[
                [ 0     0                0              0                0 ]
                [ 1/14            1/9    (13-3*√21)/63  (14- 3*√21)/126  0 ]
                [ 1/32  (91+21*√21)/576          11/72  (91-21*√21)/576  0 ]
                [ 1/14  (14+ 3*√21)/126  (13+3*√21)/63            1/9    0 ]
                [ 0               7/18            2/9             7/18   0 ]
            ]

        Tableau{T}(:LobattoIIIC̄5, 8, a, get_lobatto_weights(5), get_lobatto_nodes(5); R∞ = +1)
    end


    function _getTableauLobattoIIID2(T=Float64)
        lobC = _getTableauLobattoIIIC2(BigFloat)
        lobC̄ = _getTableauLobattoIIIC̄2(BigFloat)
        Tableau{T}(:LobattoIIID2, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c; R∞ = +1)
    end

    function _getTableauLobattoIIID3(T=Float64)
        lobC = _getTableauLobattoIIIC3(BigFloat)
        lobC̄ = _getTableauLobattoIIIC̄3(BigFloat)
        Tableau{T}(:LobattoIIID3, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c; R∞ = -1)
    end

    function _getTableauLobattoIIID4(T=Float64)
        lobC = _getTableauLobattoIIIC4(BigFloat)
        lobC̄ = _getTableauLobattoIIIC̄4(BigFloat)
        Tableau{T}(:LobattoIIID4, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c; R∞ = +1)
    end

    function _getTableauLobattoIIID5(T=Float64)
        lobC = _getTableauLobattoIIIC5(BigFloat)
        lobC̄ = _getTableauLobattoIIIC̄5(BigFloat)
        Tableau{T}(:LobattoIIID5, lobC.o, (lobC̄.a + lobC.a)/2, lobC.b, lobC.c; R∞ = -1)
    end


    function _getTableauLobattoIIIE2(T=Float64)
        lobA = _getTableauLobattoIIIA2(BigFloat)
        lobB = _getTableauLobattoIIIB2(BigFloat)
        Tableau{T}(:LobattoIIIE2, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c; R∞ = +1)
    end

    function _getTableauLobattoIIIE3(T=Float64)
        lobA = _getTableauLobattoIIIA3(BigFloat)
        lobB = _getTableauLobattoIIIB3(BigFloat)
        Tableau{T}(:LobattoIIIE3, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c; R∞ = -1)
    end

    function _getTableauLobattoIIIE4(T=Float64)
        lobA = _getTableauLobattoIIIA4(BigFloat)
        lobB = _getTableauLobattoIIIB4(BigFloat)
        Tableau{T}(:LobattoIIIE4, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c; R∞ = +1)
    end

    function _getTableauLobattoIIIE5(T=Float64)
        lobA = _getTableauLobattoIIIA5(BigFloat)
        lobB = _getTableauLobattoIIIB5(BigFloat)
        Tableau{T}(:LobattoIIIE5, lobA.o, (lobA.a + lobB.a)/2, lobA.b, lobA.c; R∞ = -1)
    end


    function _getTableauLobattoIIIF2(T=Float64)
        a = BigFloat[
                [1/12  -1/12 ]
                [7/12   5/12 ]
            ]

        Tableau{T}(:LobattoIIIF2, 4, a, get_lobatto_weights(2), get_lobatto_nodes(2); R∞ = +1)
    end

    function _getTableauLobattoIIIF3(T=Float64)
        a = BigFloat[
                [1/30   -1/15    1/30 ]
                [5/24    1/3    -1/24 ]
                [2/15   11/15    2/15 ]
            ]

        Tableau{T}(:LobattoIIIF3, 6, a, get_lobatto_weights(3), get_lobatto_nodes(3); R∞ = -1)
    end

    function _getTableauLobattoIIIF4(T=Float64)
        a = BigFloat[
                [  1/56                 -√5/56           √5/56   -1/56         ]
                [ 37/420+√5/120  5/24-   √5/210  5/24-47*√5/420  -1/210+√5/120 ]
                [ 37/420-√5/120  5/24+47*√5/420  5/24+   √5/210  -1/210-√5/120 ]
                [ 17/168         5/12-   √5/56   5/12+   √5/56   11/168        ]
            ]

        Tableau{T}(:LobattoIIIF4, 8, a, get_lobatto_weights(4), get_lobatto_nodes(4); R∞ = +1)
    end


    function _getTableauLobattoIIIG2(T=Float64)
        SymplecticTableau(_getTableauLobattoIIIF2(BigFloat); name=:LobattoIIIG2, T=T)
    end

    function _getTableauLobattoIIIG3(T=Float64)
        SymplecticTableau(_getTableauLobattoIIIF3(BigFloat); name=:LobattoIIIG3, T=T)
    end

    function _getTableauLobattoIIIG4(T=Float64)
        SymplecticTableau(_getTableauLobattoIIIF4(BigFloat); name=:LobattoIIIG4, T=T)
    end


    @test_throws ErrorException TableauLobattoIII(1)
    @test_throws ErrorException TableauLobattoIIIA(1)
    @test_throws ErrorException TableauLobattoIIIĀ(1)
    @test_throws ErrorException TableauLobattoIIIB(1)
    @test_throws ErrorException TableauLobattoIIIB̄(1)
    @test_throws ErrorException TableauLobattoIIIC(1)
    @test_throws ErrorException TableauLobattoIIIC̄(1)
    @test_throws ErrorException TableauLobattoIIID(1)
    @test_throws ErrorException TableauLobattoIIID̄(1)
    @test_throws ErrorException TableauLobattoIIIE(1)
    @test_throws ErrorException TableauLobattoIIIĒ(1)
    @test_throws ErrorException TableauLobattoIIIF(1)
    @test_throws ErrorException TableauLobattoIIIF̄(1)
    @test_throws ErrorException TableauLobattoIIIG(1)
    @test_throws ErrorException TableauLobattoIIIḠ(1)


    @test TableauLobattoIII(2) ≈ _getTableauLobattoIIIC̄2()
    @test TableauLobattoIII(3) ≈ _getTableauLobattoIIIC̄3()
    @test TableauLobattoIII(4) ≈ _getTableauLobattoIIIC̄4()
    @test TableauLobattoIII(5) ≈ _getTableauLobattoIIIC̄5()

    @test TableauLobattoIIIA(2) ≈ _getTableauLobattoIIIA2()
    @test TableauLobattoIIIA(3) ≈ _getTableauLobattoIIIA3()
    @test TableauLobattoIIIA(4) ≈ _getTableauLobattoIIIA4()
    @test TableauLobattoIIIA(5) ≈ _getTableauLobattoIIIA5()

    @test TableauLobattoIIIB(2) ≈ _getTableauLobattoIIIB2()
    @test TableauLobattoIIIB(3) ≈ _getTableauLobattoIIIB3()
    @test TableauLobattoIIIB(4) ≈ _getTableauLobattoIIIB4()
    @test TableauLobattoIIIB(5) ≈ _getTableauLobattoIIIB5()

    @test TableauLobattoIIIC(2) ≈ _getTableauLobattoIIIC2()
    @test TableauLobattoIIIC(3) ≈ _getTableauLobattoIIIC3()
    @test TableauLobattoIIIC(4) ≈ _getTableauLobattoIIIC4()
    @test TableauLobattoIIIC(5) ≈ _getTableauLobattoIIIC5()

    @test TableauLobattoIIIC̄(2) ≈ _getTableauLobattoIIIC̄2()
    @test TableauLobattoIIIC̄(3) ≈ _getTableauLobattoIIIC̄3()
    @test TableauLobattoIIIC̄(4) ≈ _getTableauLobattoIIIC̄4()
    @test TableauLobattoIIIC̄(5) ≈ _getTableauLobattoIIIC̄5()

    @test TableauLobattoIIID(2) ≈ _getTableauLobattoIIID2()
    @test TableauLobattoIIID(3) ≈ _getTableauLobattoIIID3()
    @test TableauLobattoIIID(4) ≈ _getTableauLobattoIIID4()
    @test TableauLobattoIIID(5) ≈ _getTableauLobattoIIID5()

    @test TableauLobattoIIIE(2) ≈ _getTableauLobattoIIIE2()
    @test TableauLobattoIIIE(3) ≈ _getTableauLobattoIIIE3()
    @test TableauLobattoIIIE(4) ≈ _getTableauLobattoIIIE4()
    @test TableauLobattoIIIE(5) ≈ _getTableauLobattoIIIE5()

    @test TableauLobattoIIIF(2) ≈ _getTableauLobattoIIIF2()
    @test TableauLobattoIIIF(3) ≈ _getTableauLobattoIIIF3()
    @test TableauLobattoIIIF(4) ≈ _getTableauLobattoIIIF4()

    @test TableauLobattoIIIG(2) ≈ _getTableauLobattoIIIG2()
    @test TableauLobattoIIIG(3) ≈ _getTableauLobattoIIIG3()
    @test TableauLobattoIIIG(4) ≈ _getTableauLobattoIIIG4()

    @test reference(TableauLobattoIII(2)) == reference(Val(:LobattoIII))
    @test reference(TableauLobattoIIIA(2)) == reference(Val(:LobattoIIIA))
    @test reference(TableauLobattoIIIB(2)) == reference(Val(:LobattoIIIB))
    @test reference(TableauLobattoIIIC(2)) == reference(Val(:LobattoIIIC))
    @test reference(TableauLobattoIIID(2)) == reference(Val(:LobattoIIID))
    @test reference(TableauLobattoIIIE(2)) == reference(Val(:LobattoIIIE))
    @test reference(TableauLobattoIIIF(2)) == reference(Val(:LobattoIIIF))


    for s in 2:5
        @test TableauLobattoIIIA(s) ≈ TableauLobattoIIIB̄(s)
        @test TableauLobattoIIIB(s) ≈ TableauLobattoIIIĀ(s)
        @test TableauLobattoIIID(s) ≈ TableauLobattoIIID̄(s)
        @test TableauLobattoIIIE(s) ≈ TableauLobattoIIIĒ(s)
        @test TableauLobattoIIIG(s) ≈ TableauLobattoIIIḠ(s)
    
        @test !issymplectic(TableauLobattoIII(s))
        @test !issymplectic(TableauLobattoIIIA(s))
        @test !issymplectic(TableauLobattoIIIB(s))
        @test !issymplectic(TableauLobattoIIIC(s))
        @test !issymplectic(TableauLobattoIIIC̄(s))
        @test  issymplectic(TableauLobattoIIID(s))
        @test  issymplectic(TableauLobattoIIIE(s))
    end

    for T in (Float32, Float64, BigFloat, SymPo)
        for s in 2:4
            @test_nowarn get_lobatto_nodes(T,s)
            @test_nowarn get_lobatto_weights(T,s)

            @test_nowarn get_lobatto_a_coefficients(T,s)
            @test_nowarn get_lobatto_b_coefficients(T,s)
            @test_nowarn get_lobatto_c_coefficients(T,s)
            @test_nowarn get_lobatto_c̄_coefficients(T,s)
            @test_nowarn get_lobatto_d_coefficients(T,s)
            @test_nowarn get_lobatto_e_coefficients(T,s)
            @test_nowarn get_lobatto_f_coefficients(T,s)
            @test_nowarn get_lobatto_g_coefficients(T,s)

            @test_nowarn TableauLobattoIIIA(T,s)
            @test_nowarn TableauLobattoIIIB(T,s)
            @test_nowarn TableauLobattoIIIC(T,s)
            @test_nowarn TableauLobattoIIID(T,s)
            @test_nowarn TableauLobattoIIIE(T,s)
            @test_nowarn TableauLobattoIIIF(T,s)
            @test_nowarn TableauLobattoIIIG(T,s)
        end
    end

    @test get_lobatto_nodes(Float32,2) ≈ get_lobatto_nodes(Float64,2)
    @test get_lobatto_weights(Float32,2) ≈ get_lobatto_weights(Float64,2)
    @test get_lobatto_a_coefficients(Float32,2) ≈ get_lobatto_a_coefficients(Float64,2)
    @test get_lobatto_b_coefficients(Float32,2) ≈ get_lobatto_b_coefficients(Float64,2)
    @test get_lobatto_c_coefficients(Float32,2) ≈ get_lobatto_c_coefficients(Float64,2)
    @test get_lobatto_c̄_coefficients(Float32,2) ≈ get_lobatto_c̄_coefficients(Float64,2)
    @test get_lobatto_d_coefficients(Float32,2) ≈ get_lobatto_d_coefficients(Float64,2)
    @test get_lobatto_e_coefficients(Float32,2) ≈ get_lobatto_e_coefficients(Float64,2)
    @test get_lobatto_f_coefficients(Float32,2) ≈ get_lobatto_f_coefficients(Float64,2)
    @test get_lobatto_g_coefficients(Float32,2) ≈ get_lobatto_g_coefficients(Float64,2)

    @test get_lobatto_nodes(SymPo,2) ≈ get_lobatto_nodes(Float64,2)
    @test get_lobatto_weights(SymPo,2) ≈ get_lobatto_weights(Float64,2)
    @test get_lobatto_a_coefficients(SymPo,2) ≈ get_lobatto_a_coefficients(Float64,2)
    @test get_lobatto_b_coefficients(SymPo,2) ≈ get_lobatto_b_coefficients(Float64,2)
    @test get_lobatto_c_coefficients(SymPo,2) ≈ get_lobatto_c_coefficients(Float64,2)
    @test get_lobatto_c̄_coefficients(SymPo,2) ≈ get_lobatto_c̄_coefficients(Float64,2)
    @test get_lobatto_d_coefficients(SymPo,2) ≈ get_lobatto_d_coefficients(Float64,2)
    @test get_lobatto_e_coefficients(SymPo,2) ≈ get_lobatto_e_coefficients(Float64,2)
    @test get_lobatto_f_coefficients(SymPo,2) ≈ get_lobatto_f_coefficients(Float64,2)
    @test get_lobatto_g_coefficients(SymPo,2) ≈ get_lobatto_g_coefficients(Float64,2)

    @test TableauLobattoIIIA(Float32,2) ≈ TableauLobattoIIIA(Float64,2)
    @test TableauLobattoIIIB(Float32,2) ≈ TableauLobattoIIIB(Float64,2)
    @test TableauLobattoIIIC(Float32,2) ≈ TableauLobattoIIIC(Float64,2)
    @test TableauLobattoIIID(Float32,2) ≈ TableauLobattoIIID(Float64,2)
    @test TableauLobattoIIIE(Float32,2) ≈ TableauLobattoIIIE(Float64,2)
    @test TableauLobattoIIIF(Float32,2) ≈ TableauLobattoIIIF(Float64,2)
    @test TableauLobattoIIIG(Float32,2) ≈ TableauLobattoIIIG(Float64,2)
    
    @test TableauLobattoIIIA(SymPo,2) ≈ TableauLobattoIIIA(Float64,2)
    @test TableauLobattoIIIB(SymPo,2) ≈ TableauLobattoIIIB(Float64,2)
    @test TableauLobattoIIIC(SymPo,2) ≈ TableauLobattoIIIC(Float64,2)
    @test TableauLobattoIIID(SymPo,2) ≈ TableauLobattoIIID(Float64,2)
    @test TableauLobattoIIIE(SymPo,2) ≈ TableauLobattoIIIE(Float64,2)
    @test TableauLobattoIIIF(SymPo,2) ≈ TableauLobattoIIIF(Float64,2)
    @test TableauLobattoIIIG(SymPo,2) ≈ TableauLobattoIIIG(Float64,2)

end
