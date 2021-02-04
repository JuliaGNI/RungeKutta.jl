using RungeKutta: name, order, nstages, coefficients, weights, nodes

@testset "$(rpad("Partitioned Tableaus",80))" begin

    @test typeof(PartitionedTableau(:PRK4, TableauRK4(), TableauRK4())) <: PartitionedTableau
    @test PartitionedTableau(:PRK4, TableauRK4(), TableauRK4()) == PartitionedTableau(:PRK4, TableauRK4())
    

    @test typeof(PartitionedTableauGauss(1)) <: PartitionedTableau
    @test order(PartitionedTableauGauss(1)) == 2
    @test nstages(PartitionedTableauGauss(1)) == 1
    @test issymplectic(PartitionedTableauGauss(1))
    @test issymmetric(PartitionedTableauGauss(1))
    @test all(all.(check_symmetry(PartitionedTableauGauss(1))))

    @test !isexplicit(PartitionedTableauGauss(1))
    @test  isimplicit(PartitionedTableauGauss(1))
    @test !isdiagnonallyimplicit(PartitionedTableauGauss(1))
    @test  isfullyimplicit(PartitionedTableauGauss(1))
    
    @test typeof(PartitionedTableauGauss(2)) <: PartitionedTableau
    @test order(PartitionedTableauGauss(2)) == 4
    @test nstages(PartitionedTableauGauss(2)) == 2
    @test issymplectic(PartitionedTableauGauss(2))
    @test issymmetric(PartitionedTableauGauss(2))
    @test all(all.(check_symmetry(PartitionedTableauGauss(2))))

    @test !isexplicit(PartitionedTableauGauss(2))
    @test  isimplicit(PartitionedTableauGauss(2))
    @test !isdiagnonallyimplicit(PartitionedTableauGauss(2))
    @test  isfullyimplicit(PartitionedTableauGauss(2))
    
    
    @test typeof(TableauLobattoIIIAIIIB(2)) <: PartitionedTableau
    @test order(TableauLobattoIIIAIIIB(2)) == 2
    @test nstages(TableauLobattoIIIAIIIB(2)) == 2
    @test issymplectic(TableauLobattoIIIAIIIB(2))
    @test issymmetric(TableauLobattoIIIAIIIB(2))

    @test  isexplicit(TableauLobattoIIIAIIIB(2))
    @test !isimplicit(TableauLobattoIIIAIIIB(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIIAIIIB(2))
    @test !isfullyimplicit(TableauLobattoIIIAIIIB(2))
    
    @test typeof(TableauLobattoIIIAIIIB(3)) <: PartitionedTableau
    @test order(TableauLobattoIIIAIIIB(3)) == 4
    @test nstages(TableauLobattoIIIAIIIB(3)) == 3
    @test issymplectic(TableauLobattoIIIAIIIB(3))
    @test issymmetric(TableauLobattoIIIAIIIB(3))

    @test !isexplicit(TableauLobattoIIIAIIIB(3))
    @test  isimplicit(TableauLobattoIIIAIIIB(3))
    @test !isdiagnonallyimplicit(TableauLobattoIIIAIIIB(3))
    @test  isfullyimplicit(TableauLobattoIIIAIIIB(3))
    

    @test typeof(TableauLobattoIIIBIIIA(2)) <: PartitionedTableau
    @test order(TableauLobattoIIIBIIIA(2)) == 2
    @test nstages(TableauLobattoIIIBIIIA(2)) == 2
    @test issymplectic(TableauLobattoIIIBIIIA(2))
    @test issymmetric(TableauLobattoIIIBIIIA(2))

    @test  isexplicit(TableauLobattoIIIBIIIA(2))
    @test !isimplicit(TableauLobattoIIIBIIIA(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIIBIIIA(2))
    @test !isfullyimplicit(TableauLobattoIIIBIIIA(2))
    
    @test typeof(TableauLobattoIIIBIIIA(3)) <: PartitionedTableau
    @test order(TableauLobattoIIIBIIIA(3)) == 4
    @test nstages(TableauLobattoIIIBIIIA(3)) == 3
    @test issymplectic(TableauLobattoIIIBIIIA(3))
    @test issymmetric(TableauLobattoIIIBIIIA(3))

    @test !isexplicit(TableauLobattoIIIBIIIA(3))
    @test  isimplicit(TableauLobattoIIIBIIIA(3))
    @test !isdiagnonallyimplicit(TableauLobattoIIIBIIIA(3))
    @test  isfullyimplicit(TableauLobattoIIIBIIIA(3))
    

    @test typeof(TableauLobattoIIICIIIC̄(2)) <: PartitionedTableau
    @test order(TableauLobattoIIICIIIC̄(2)) == 2
    @test nstages(TableauLobattoIIICIIIC̄(2)) == 2
    @test  issymplectic(TableauLobattoIIICIIIC̄(2))
    @test !issymmetric(TableauLobattoIIICIIIC̄(2))

    @test !isexplicit(TableauLobattoIIICIIIC̄(2))
    @test  isimplicit(TableauLobattoIIICIIIC̄(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIICIIIC̄(2))
    @test  isfullyimplicit(TableauLobattoIIICIIIC̄(2))
    
    @test typeof(TableauLobattoIIICIIIC̄(3)) <: PartitionedTableau
    @test order(TableauLobattoIIICIIIC̄(3)) == 4
    @test nstages(TableauLobattoIIICIIIC̄(3)) == 3
    @test  issymplectic(TableauLobattoIIICIIIC̄(3))
    @test !issymmetric(TableauLobattoIIICIIIC̄(3))

    @test !isexplicit(TableauLobattoIIICIIIC̄(3))
    @test  isimplicit(TableauLobattoIIICIIIC̄(3))
    @test !isdiagnonallyimplicit(TableauLobattoIIICIIIC̄(3))
    @test  isfullyimplicit(TableauLobattoIIICIIIC̄(3))
    
    
    @test typeof(TableauLobattoIIIC̄IIIC(2)) <: PartitionedTableau
    @test order(TableauLobattoIIIC̄IIIC(2)) == 2
    @test nstages(TableauLobattoIIIC̄IIIC(2)) == 2
    @test  issymplectic(TableauLobattoIIIC̄IIIC(2))
    @test !issymmetric(TableauLobattoIIIC̄IIIC(2))

    @test !isexplicit(TableauLobattoIIIC̄IIIC(2))
    @test  isimplicit(TableauLobattoIIIC̄IIIC(2))
    @test !isdiagnonallyimplicit(TableauLobattoIIIC̄IIIC(2))
    @test  isfullyimplicit(TableauLobattoIIIC̄IIIC(2))
    
    @test typeof(TableauLobattoIIIC̄IIIC(3)) <: PartitionedTableau
    @test order(TableauLobattoIIIC̄IIIC(3)) == 4
    @test nstages(TableauLobattoIIIC̄IIIC(3)) == 3
    @test  issymplectic(TableauLobattoIIIC̄IIIC(3))
    @test !issymmetric(TableauLobattoIIIC̄IIIC(3))

    @test !isexplicit(TableauLobattoIIIC̄IIIC(3))
    @test  isimplicit(TableauLobattoIIIC̄IIIC(3))
    @test !isdiagnonallyimplicit(TableauLobattoIIIC̄IIIC(3))
    @test  isfullyimplicit(TableauLobattoIIIC̄IIIC(3))

end
