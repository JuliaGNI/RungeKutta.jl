using RungeKutta
using Test

@testset "$(rpad("Symmetry",80))" begin
    for s in 1:10
        @test check_symmetry(TableauGauss(s)) == Array{Bool}(ones(s, s))
        @test issymmetric(TableauGauss(s))
    end

    for s in 2:10
        @test check_symmetry(TableauLobattoIIIA(s)) == Array{Bool}(ones(s, s))
        @test check_symmetry(TableauLobattoIIIB(s)) == Array{Bool}(ones(s, s))
        @test check_symmetry(TableauLobattoIIIC(s)) != Array{Bool}(ones(s, s))
        @test check_symmetry(TableauLobattoIIIC̄(s)) != Array{Bool}(ones(s, s))
        @test check_symmetry(TableauLobattoIIID(s)) == Array{Bool}(ones(s, s))
        @test check_symmetry(TableauLobattoIIIE(s)) == Array{Bool}(ones(s, s))
        @test check_symmetry(TableauRadauIA(s)) != Array{Bool}(ones(s, s))
        @test check_symmetry(TableauRadauIIA(s)) != Array{Bool}(ones(s, s))

        @test issymmetric(TableauLobattoIIIA(s))
        @test issymmetric(TableauLobattoIIIB(s))
        @test !issymmetric(TableauLobattoIIIC(s))
        @test !issymmetric(TableauLobattoIIIC̄(s))
        @test issymmetric(TableauLobattoIIID(s))
        @test issymmetric(TableauLobattoIIIE(s))
        @test !issymmetric(TableauRadauIA(s))
        @test !issymmetric(TableauRadauIIA(s))
    end
end

@testset "$(rpad("Symmetry in Rational{BigInt}",80))" begin
    QQ = Rational{BigInt}

    # the default tolerance of an exact element type is exactly zero
    midpoint = Tableau(:midpoint, 2, QQ[1//2;;], QQ[1], QQ[1 // 2])
    trapezoidal = Tableau(:trapezoidal, 2, QQ[0 0; 1//2 1//2], QQ[1 // 2, 1 // 2], QQ[0, 1])
    @test eltype(midpoint) == eltype(trapezoidal) == QQ
    @test check_symmetry(midpoint) == trues(1, 1)
    @test check_symmetry(trapezoidal) == trues(2, 2)
    @test issymmetric(midpoint)
    @test issymmetric(trapezoidal)
    @test !issymmetric(TableauRK4(QQ))
    @test check_symmetry(PartitionedTableau(:trapezoidal, trapezoidal)) ==
          (trues(2, 2), trues(2, 2))
end
