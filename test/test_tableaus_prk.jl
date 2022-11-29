using RungeKutta: name, order, nstages, coefficients, weights, nodes

@testset "$(rpad("Partitioned Tableaus",80))" begin

    @test typeof(PartitionedTableau(:PRK4, TableauRK4(), TableauRK4())) <: PartitionedTableau
    @test PartitionedTableau(:PRK4, TableauRK4(), TableauRK4()) == PartitionedTableau(:PRK4, TableauRK4())
    

    @test typeof(PartitionedTableauGauss(1)) <: PartitionedTableau
    @test typeof(PartitionedTableauGauss(2)) <: PartitionedTableau

    @test PartitionedTableauGauss(1) == PartitionedTableauGauss(Float64, 1)
    @test PartitionedTableauGauss(2) == PartitionedTableauGauss(Float64, 2)

    @test order(PartitionedTableauGauss(1)) == 2
    @test order(PartitionedTableauGauss(2)) == 4

    @test nstages(PartitionedTableauGauss(1)) == 1
    @test nstages(PartitionedTableauGauss(2)) == 2

    @test issymplectic(PartitionedTableauGauss(1))
    @test issymplectic(PartitionedTableauGauss(2))

    @test issymmetric(PartitionedTableauGauss(1))
    @test issymmetric(PartitionedTableauGauss(2))

    @test all(all.(check_symmetry(PartitionedTableauGauss(1))))
    @test all(all.(check_symmetry(PartitionedTableauGauss(2))))

    @test PartitionedTableauGauss(1).R∞ == -1
    @test PartitionedTableauGauss(2).R∞ == +1

    @test !isexplicit(PartitionedTableauGauss(1))
    @test !isexplicit(PartitionedTableauGauss(2))

    @test  isimplicit(PartitionedTableauGauss(1))
    @test  isimplicit(PartitionedTableauGauss(2))

    @test !isdiagnonallyimplicit(PartitionedTableauGauss(1))
    @test !isdiagnonallyimplicit(PartitionedTableauGauss(2))

    @test  isfullyimplicit(PartitionedTableauGauss(1))
    @test  isfullyimplicit(PartitionedTableauGauss(2))
    
    
    @test typeof(TableauLobattoIIIAIIIB(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIAIIIB(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIBIIIA(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIBIIIA(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIAIIIĀ(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIAIIIĀ(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIBIIIB̄(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIBIIIB̄(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIICIIIC̄(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIICIIIC̄(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIC̄IIIC(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIC̄IIIC(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIDIIID̄(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIDIIID̄(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIEIIIĒ(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIEIIIĒ(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIFIIIF̄(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIFIIIF̄(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIF̄IIIF(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIF̄IIIF(3)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIGIIIḠ(2)) <: PartitionedTableau
    @test typeof(TableauLobattoIIIGIIIḠ(3)) <: PartitionedTableau

    @test TableauLobattoIIIAIIIB(2) == TableauLobattoIIIAIIIB(Float64, 2)
    @test TableauLobattoIIIBIIIA(2) == TableauLobattoIIIBIIIA(Float64, 2)
    @test TableauLobattoIIIAIIIĀ(2) == TableauLobattoIIIAIIIĀ(Float64, 2)
    @test TableauLobattoIIIBIIIB̄(2) == TableauLobattoIIIBIIIB̄(Float64, 2)
    @test TableauLobattoIIICIIIC̄(2) == TableauLobattoIIICIIIC̄(Float64, 2)
    @test TableauLobattoIIIC̄IIIC(2) == TableauLobattoIIIC̄IIIC(Float64, 2)
    @test TableauLobattoIIIDIIID̄(2) == TableauLobattoIIIDIIID̄(Float64, 2)
    @test TableauLobattoIIIEIIIĒ(2) == TableauLobattoIIIEIIIĒ(Float64, 2)
    @test TableauLobattoIIIFIIIF̄(2) == TableauLobattoIIIFIIIF̄(Float64, 2)
    @test TableauLobattoIIIF̄IIIF(2) == TableauLobattoIIIF̄IIIF(Float64, 2)
    @test TableauLobattoIIIGIIIḠ(2) == TableauLobattoIIIGIIIḠ(Float64, 2)

    @test order(TableauLobattoIIIAIIIB(2)) == 2
    @test order(TableauLobattoIIIAIIIB(3)) == 4
    @test order(TableauLobattoIIIBIIIA(2)) == 2
    @test order(TableauLobattoIIIBIIIA(3)) == 4
    @test order(TableauLobattoIIICIIIC̄(2)) == 2
    @test order(TableauLobattoIIICIIIC̄(3)) == 4
    @test order(TableauLobattoIIIC̄IIIC(2)) == 2
    @test order(TableauLobattoIIIC̄IIIC(3)) == 4

    @test nstages(TableauLobattoIIIAIIIB(2)) == 2
    @test nstages(TableauLobattoIIIAIIIB(3)) == 3
    @test nstages(TableauLobattoIIIBIIIA(2)) == 2
    @test nstages(TableauLobattoIIIBIIIA(3)) == 3
    @test nstages(TableauLobattoIIICIIIC̄(2)) == 2
    @test nstages(TableauLobattoIIICIIIC̄(3)) == 3
    @test nstages(TableauLobattoIIIC̄IIIC(2)) == 2
    @test nstages(TableauLobattoIIIC̄IIIC(3)) == 3

    @test issymplectic(TableauLobattoIIIAIIIB(2))
    @test issymplectic(TableauLobattoIIIAIIIB(3))
    @test issymplectic(TableauLobattoIIIBIIIA(2))
    @test issymplectic(TableauLobattoIIIBIIIA(3))
    @test issymplectic(TableauLobattoIIICIIIC̄(2))
    @test issymplectic(TableauLobattoIIICIIIC̄(3))
    @test issymplectic(TableauLobattoIIIC̄IIIC(2))
    @test  issymplectic(TableauLobattoIIIC̄IIIC(3))

    @test  issymmetric(TableauLobattoIIIAIIIB(2))
    @test  issymmetric(TableauLobattoIIIAIIIB(3))
    @test  issymmetric(TableauLobattoIIIBIIIA(2))
    @test  issymmetric(TableauLobattoIIIBIIIA(3))
    @test !issymmetric(TableauLobattoIIICIIIC̄(2))
    @test !issymmetric(TableauLobattoIIICIIIC̄(3))
    @test !issymmetric(TableauLobattoIIIC̄IIIC(2))
    @test !issymmetric(TableauLobattoIIIC̄IIIC(3))

    @test TableauLobattoIIIAIIIB(2).R∞ == -1
    @test TableauLobattoIIIAIIIB(3).R∞ == +1
    @test TableauLobattoIIIBIIIA(2).R∞ == -1
    @test TableauLobattoIIIBIIIA(3).R∞ == +1
    @test TableauLobattoIIICIIIC̄(2).R∞ == -1
    @test TableauLobattoIIICIIIC̄(3).R∞ == +1
    @test TableauLobattoIIIC̄IIIC(2).R∞ == -1
    @test TableauLobattoIIIC̄IIIC(3).R∞ == +1

    @test  isexplicit(TableauLobattoIIIAIIIB(2))
    @test !isexplicit(TableauLobattoIIIAIIIB(3))
    @test  isexplicit(TableauLobattoIIIBIIIA(2))
    @test !isexplicit(TableauLobattoIIIBIIIA(3))
    @test !isexplicit(TableauLobattoIIICIIIC̄(2))
    @test !isexplicit(TableauLobattoIIICIIIC̄(3))
    @test !isexplicit(TableauLobattoIIIC̄IIIC(2))
    @test !isexplicit(TableauLobattoIIIC̄IIIC(3))

    @test !isimplicit(TableauLobattoIIIAIIIB(2))
    @test  isimplicit(TableauLobattoIIIAIIIB(3))
    @test !isimplicit(TableauLobattoIIIBIIIA(2))
    @test  isimplicit(TableauLobattoIIIBIIIA(3))
    @test  isimplicit(TableauLobattoIIICIIIC̄(2))
    @test  isimplicit(TableauLobattoIIICIIIC̄(3))
    @test  isimplicit(TableauLobattoIIIC̄IIIC(2))
    @test  isimplicit(TableauLobattoIIIC̄IIIC(3))

    @test !isdiagnonallyimplicit(TableauLobattoIIIAIIIB(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIIAIIIB(3))
    @test !isdiagnonallyimplicit(TableauLobattoIIIBIIIA(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIIBIIIA(3))
    @test !isdiagnonallyimplicit(TableauLobattoIIICIIIC̄(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIICIIIC̄(3))
    @test !isdiagnonallyimplicit(TableauLobattoIIIC̄IIIC(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIIC̄IIIC(3))

    @test !isfullyimplicit(TableauLobattoIIIAIIIB(2))
    @test  isfullyimplicit(TableauLobattoIIIAIIIB(3))
    @test !isfullyimplicit(TableauLobattoIIIBIIIA(2))
    @test  isfullyimplicit(TableauLobattoIIIBIIIA(3))
    @test  isfullyimplicit(TableauLobattoIIICIIIC̄(2))
    @test  isfullyimplicit(TableauLobattoIIICIIIC̄(3))
    @test  isfullyimplicit(TableauLobattoIIIC̄IIIC(2))
    @test  isfullyimplicit(TableauLobattoIIIC̄IIIC(3))

end
