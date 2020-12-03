import LinearAlgebra: I
import RungeKutta: @big, _big, istriustrict, istrilstrict

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


    for n in 1:5
        A = zeros(n,n)
        B = zeros(n,n)

        for i in 1:n
            for j in i+1:n
                A[i,j] = rand()
                B[j,i] = rand()
            end
        end

        @test istriustrict(A)
        @test istrilstrict(B)

        @test !istriustrict(A .+ Matrix(I, n, n))
        @test !istrilstrict(B .+ Matrix(I, n, n))

        @test !istriustrict(rand(n,n))
        @test !istrilstrict(rand(n,n))
    end

end
