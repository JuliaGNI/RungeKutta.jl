import LinearAlgebra: I, norm, nullspace
import Polynomials: Polynomial
import RungeKutta: istriustrict, istrilstrict
import RungeKutta: _legendre, _shifted_legendre, _nullvector

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


    for T in (Float32, Float64, BigFloat, symtype())
        for s in 1:3
            @test_nowarn _legendre(s, Polynomial(T[0,1]))
            @test_nowarn _shifted_legendre(s,T)
        end
    end


    # The nullvector of a rank-deficient matrix agrees with the singular value
    # decomposition, and is computed generically, i.e. also for BigFloat.
    for T in (Float64, BigFloat)
        A = T[1 2 3; 2 4 6; 1 1 1]     # rank 2, nullvector ∝ [1, -2, 1]
        w = _nullvector(A)

        @test eltype(w) == T
        @test norm(w) ≈ one(T)
        @test norm(A * w) < 16eps(T)
        @test abs.(w) ≈ abs.(nullspace(Float64.(A))[:,begin])
    end

end
