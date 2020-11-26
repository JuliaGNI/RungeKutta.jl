@testset "$(rpad("Radau Tableaus",80))" begin

    function _TableauRadauIIA2(T=Float64)
        a = [[5//12  -1//12]
             [3//4    1//4 ]]
        b = [3//4, 1//4]
        c = [1//3, 1//1]

        Tableau{T}(:RadauIIA2, 3, a, b, c)
    end

    function _TableauRadauIIA3(T=Float64)
        a = [
                [11/45 -  7*√6/360    37/225-169*√6/1800   -2/225+√6/75]
                [37/225+169*√6/1800   11/45 +  7*√6/360    -2/225-√6/75]
                [ 4/9  -    √6/36      4/9  +    √6/36      1/9        ]
            ]
        b = [4/9-√6/36,   4/9+√6/36,  1/9 ]
        c = [2/5-√6/10,   2/5+√6/10,  1   ]

        Tableau{T}(:RadauIIA3, 5, a, b, c)
    end

    @test TableauRadauIIA(2) ≈ _TableauRadauIIA2()
    @test TableauRadauIIA(3) ≈ _TableauRadauIIA3()

end
