@testset "$(rpad("Symplecticity",80))" begin

    for s in 1:10
        g = TableauGauss(s)
        g̃ = SymplecticConjugateTableau(g)

        @test g.a ≈ g̃.a atol=1E-15
        @test g.b == g̃.b
        @test g.c == g̃.c

        @test symplecticity_error(g) ≈ zeros(s,s)  atol=eps()
        @test check_symplecticity(g) == Array{Bool}(ones(s,s))

        @test issymplectic(g)
        @test issymplectic(g̃)
        @test issymplectic(SymplecticTableau(g))
        @test issymplectic(SymplecticPartitionedTableau(g))
    end

    for s in 2:10
        A = TableauLobattoIIIA(s)
        B = TableauLobattoIIIB(s)
        Ã = SymplecticConjugateTableau(A)
        B̃ = SymplecticConjugateTableau(B)

        @test A.a ≈ B̃.a atol=1E-15
        @test A.b == Ã.b
        @test A.c == Ã.c

        @test B.a ≈ Ã.a atol=1E-15
        @test B.b == Ã.b
        @test B.c == Ã.c


        Â = SymplecticTableau(A)
        B̂ = SymplecticTableau(B)
        E = TableauLobattoIIIE(s)

        @test Â.a ≈ E.a atol=1E-15
        @test Â.b == E.b
        @test Â.c == E.c

        @test B̂.a ≈ E.a atol=1E-15
        @test B̂.b == E.b
        @test B̂.c == E.c

        @test issymplectic(Â)
        @test issymplectic(B̂)
        @test issymplectic(E)
        @test issymplectic(SymplecticPartitionedTableau(A))
        @test issymplectic(SymplecticPartitionedTableau(B))
    end

end
