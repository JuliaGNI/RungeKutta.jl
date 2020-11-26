@testset "$(rpad("Symplecticity",80))" begin

    for s in 1:5
        g = TableauGauss(s)
        g̃ = get_symplectic_conjugate_coefficients(g)

        @test g.a ≈ g̃.a atol=1E-15
        @test g.b == g̃.b
        @test g.c == g̃.c

        @test compute_symplecticity_error(g) ≈ zeros(s,s)  atol=eps()
        @test check_symplecticity(g) == Array{Bool}(ones(s,s))
        @test issymplectic(g)
    end

    for s in 2:5
        A = TableauLobattoIIIA(s)
        B = TableauLobattoIIIB(s)
        E = TableauLobattoIIIE(s)
        Ã = get_symplectic_conjugate_coefficients(A)
        B̃ = get_symplectic_conjugate_coefficients(B)

        @test A.a ≈ B̃.a atol=1E-15
        @test A.b == Ã.b
        @test A.c == Ã.c

        @test B.a ≈ Ã.a atol=1E-15
        @test B.b == Ã.b
        @test B.c == Ã.c

        Â = symplecticize(A)
        B̂ = symplecticize(B)

        @test Â.a ≈ E.a atol=1E-15
        @test Â.b == E.b
        @test Â.c == E.c

        @test B̂.a ≈ E.a atol=1E-15
        @test B̂.b == E.b
        @test B̂.c == E.c

        @test !issymplectic(A)
        @test !issymplectic(B)
        @test issymplectic(Â)
        @test issymplectic(B̂)
        @test issymplectic(E)
    end

end
