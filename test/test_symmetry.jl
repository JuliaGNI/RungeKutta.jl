@testset "$(rpad("Symmetry",80))" begin

    for s in 1:5
        @test check_symmetry(TableauGauss(s)) == Array{Bool}(ones(s,s))
    end

    for s in 6:10
        @test check_symmetry(TableauGauss(s); atol=1E-11, rtol=1E-11) == Array{Bool}(ones(s,s))
    end

    for s in 2:10
        @test check_symmetry(TableauLobattoIIIA(s)) == Array{Bool}(ones(s,s))
        @test check_symmetry(TableauLobattoIIIB(s)) == Array{Bool}(ones(s,s))
        @test check_symmetry(TableauLobattoIIIC(s)) != Array{Bool}(ones(s,s))
        @test check_symmetry(TableauLobattoIIIC̄(s)) != Array{Bool}(ones(s,s))
        @test check_symmetry(TableauLobattoIIID(s)) == Array{Bool}(ones(s,s))
        @test check_symmetry(TableauLobattoIIIE(s)) == Array{Bool}(ones(s,s))
    end

end
