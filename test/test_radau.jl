using RungeKutta.Tableaus: get_radau_1_nodes, get_radau_1_weights, get_radau_1_coefficients,
                           get_radau_2_nodes, get_radau_2_weights, get_radau_2_coefficients

@testset "$(rpad("Radau Tableaus",80))" begin

    @test_throws ErrorException get_radau_1_nodes(1)
    @test_throws ErrorException get_radau_1_weights(1)
    @test_throws ErrorException get_radau_1_coefficients(1)

    @test_throws ErrorException get_radau_2_nodes(1)
    @test_throws ErrorException get_radau_2_weights(1)
    @test_throws ErrorException get_radau_2_coefficients(1)

    
    function _TableauRadauIA2(T=Float64)
        a = [[1//4  -1//4  ]
             [1//4   5//12 ]]
        b = [1//4, 3//4]
        c = [0,    2//3]

        Tableau{T}(:RadauIA2, 3, a, b, c)
    end

    function _TableauRadauIA3(T=Float64)
        a = [
                [ 1//9     (- 1 -    √6)/18    (- 1 +    √6)/18  ]
                [ 1//9     ( 88 +  7*√6)/360   ( 88 - 43*√6)/360 ]
                [ 1//9     ( 88 + 43*√6)/360   ( 88 -  7*√6)/360 ]
            ]
        b = [1/9,         (16+√6)/36,  (16-√6)/36 ]
        c = [0,           ( 6-√6)/10,  ( 6+√6)/10 ]

        Tableau{T}(:RadauIA3, 5, a, b, c)
    end

    function _TableauRadauIB2(T=Float64)
        a = [[1//8  -1//8  ]
             [7//24  3//8  ]]
        b = [1//4, 3//4]
        c = [0,    2//3]

        Tableau{T}(:RadauIB2, 3, a, b, c)
    end

    function _TableauRadauIB3(T=Float64)
        a = [
                [            1//18    (- 1 -     √6)/36    (- 1 +     √6)/36   ]
                [ ( 52 +  3*√6)/450   ( 16 +     √6)/72    (472 - 217*√6)/1800 ]
                [ ( 52 -  3*√6)/450   (472 + 217*√6)/1800  ( 16 -     √6)/72   ]
            ]
        b = [1/9,         (16+√6)/36,  (16-√6)/36 ]
        c = [0,           ( 6-√6)/10,  ( 6+√6)/10 ]

        Tableau{T}(:RadauIB3, 5, a, b, c)
    end

    function _TableauRadauIIA2(T=Float64)
        a = [[5//12  -1//12]
             [3//4    1//4 ]]
        b = [3//4, 1//4]
        c = [1//3, 1//1]

        Tableau{T}(:RadauIIA2, 3, a, b, c)
    end

    function _TableauRadauIIA3(T=Float64)
        a = [
                [ 11/45 -  7*√6/360    37/225-169*√6/1800   -2/225+√6/75 ]
                [ 37/225+169*√6/1800   11/45 +  7*√6/360    -2/225-√6/75 ]
                [  4/9  -    √6/36      4/9  +    √6/36      1/9         ]
            ]
        b = [4/9-√6/36,   4/9+√6/36,  1/9 ]
        c = [2/5-√6/10,   2/5+√6/10,  1   ]

        Tableau{T}(:RadauIIA3, 5, a, b, c)
    end

    function _TableauRadauIIB2(T=Float64)
        a = [[3//8   -1//24]
             [7//8    1//8 ]]
        b = [3//4, 1//4]
        c = [1//3, 1//1]

        Tableau{T}(:RadauIIB2, 3, a, b, c)
    end

    function _TableauRadauIIB3(T=Float64)
        a = [
                [( 16 -     √6)/72    (328 - 167*√6)/1800   (-2 + 3*√6)/450 ]
                [(328 + 167*√6)/1800  ( 16 +     √6)/72     (-2 - 3*√6)/450 ]
                [( 85 -  10*√6)/180   ( 85 +  10*√6)/180              1/18  ]
            ]
        b = [4/9-√6/36,   4/9+√6/36,  1/9 ]
        c = [2/5-√6/10,   2/5+√6/10,  1   ]

        Tableau{T}(:RadauIIB3, 5, a, b, c)
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

end
