import LinearAlgebra: I, norm, nullspace
import RungeKutta: istriustrict, istrilstrict, _nullvector

@testset "$(rpad("Utility Functions",80))" begin
    for n in 1:5
        A = zeros(n, n)
        B = zeros(n, n)

        for i in 1:n
            for j in (i + 1):n
                A[i, j] = rand()
                B[j, i] = rand()
            end
        end

        @test istriustrict(A)
        @test istrilstrict(B)

        @test !istriustrict(A .+ Matrix(I, n, n))
        @test !istrilstrict(B .+ Matrix(I, n, n))

        @test !istriustrict(rand(n, n))
        @test !istrilstrict(rand(n, n))
    end

    # The nullvector of a rank-deficient matrix agrees with the singular value
    # decomposition, and is computed generically, i.e. also for BigFloat.
    for T in (Float64, BigFloat)
        A = T[1 2 3; 2 4 6; 1 1 1]     # rank 2, nullvector ∝ [1, -2, 1]
        w = _nullvector(A)

        @test eltype(w) == T
        @test norm(w) ≈ one(T)
        @test norm(A * w) < 16eps(T)
        @test abs.(w) ≈ abs.(nullspace(Float64.(A))[:, begin])

        # The sign is canonical, so the vector does not depend on the pivot order:
        # permuting the columns of A permutes w and leaves the sign alone.
        @test w[begin] > 0
        @test _nullvector(A[:, [3, 2, 1]]) ≈ reverse(w)
    end

    # A matrix whose nullspace is not one-dimensional is rejected rather than
    # answered with a vector that does not span anything.
    @test_throws ArgumentError _nullvector(Float64[1 2 3; 4 5 6; 7 8 10])   # rank 3
    @test_throws ArgumentError _nullvector(zeros(3, 3))                      # rank 0
end
