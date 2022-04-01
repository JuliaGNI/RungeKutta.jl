import LinearAlgebra: I
import RungeKutta: istriustrict, istrilstrict

@testset "$(rpad("Utility Functions",80))" begin

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
