import RungeKutta: @big, _big

@testset "$(rpad("Utility Functions",80))" begin
    
    @test _big(1)    == BigFloat(1)
    @test _big(1//1) == BigFloat(1)
    @test _big(1.0)  == BigFloat(1)
    @test _big("1")  == BigFloat(1)
    @test _big(:(1)) == BigFloat(1)

    x = @big [1  2.0  1//3]
    y = big.([1  2.0  1//3])
    z = [big(1)  big(2.0)  big(1//3)]

    @test x != y
    @test x == z

end
