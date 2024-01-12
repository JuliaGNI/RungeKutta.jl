import LinearAlgebra: I
import Polynomials: Polynomial
import RungeKutta: istriustrict, istrilstrict
import RungeKutta: _legendre, _shifted_legendre

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


    for T in (Float32, Float64, BigFloat, SymP)
        for s in 1:3
            @test_nowarn _legendre(s, Polynomial(T[0,1]))
            @test_nowarn _shifted_legendre(s,T)
        end
    end
    
end
