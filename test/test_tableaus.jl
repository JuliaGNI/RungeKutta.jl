using RungeKutta: name, order, nstages, coefficients, weights, nodes

@testset "$(rpad("ERK Tableaus",80))" begin

    @test typeof(TableauExplicitEuler()) <: Tableau
    @test order(TableauExplicitEuler()) == 1
    @test nstages(TableauExplicitEuler()) == 1

    @test  isexplicit(TableauExplicitEuler())
    @test !isimplicit(TableauExplicitEuler())
    @test !isdiagnonallyimplicit(TableauExplicitEuler())
    @test !isfullyimplicit(TableauExplicitEuler())
    
    @test TableauExplicitEuler() == TableauExplicitEuler()

    @test typeof(TableauExplicitMidpoint()) <: Tableau
    @test order(TableauExplicitMidpoint()) == 2
    @test nstages(TableauExplicitMidpoint()) == 2

    @test  isexplicit(TableauExplicitMidpoint())
    @test !isimplicit(TableauExplicitMidpoint())
    @test !isdiagnonallyimplicit(TableauExplicitMidpoint())
    @test !isfullyimplicit(TableauExplicitMidpoint())

    @test typeof(TableauHeun2()) <: Tableau
    @test order(TableauHeun2()) == 2
    @test nstages(TableauHeun2()) == 2

    @test typeof(TableauHeun3()) <: Tableau
    @test order(TableauHeun3()) == 3
    @test nstages(TableauHeun3()) == 3

    @test typeof(TableauRalston2()) <: Tableau
    @test order(TableauRalston2()) == 2
    @test nstages(TableauRalston2()) == 2

    @test typeof(TableauRalston3()) <: Tableau
    @test order(TableauRalston3()) == 3
    @test nstages(TableauRalston3()) == 3

    @test typeof(TableauRunge()) <: Tableau
    @test order(TableauRunge()) == 2
    @test nstages(TableauRunge()) == 2

    @test TableauRunge2() == TableauRunge()

    @test typeof(TableauKutta()) <: Tableau
    @test order(TableauKutta()) == 3
    @test nstages(TableauKutta()) == 3

    @test TableauKutta3() == TableauKutta()

    @test typeof(TableauRK416()) <: Tableau
    @test order(TableauRK416()) == 4
    @test nstages(TableauRK416()) == 4

    @test TableauRK4() == TableauRK416()

    @test  isexplicit(TableauRK4())
    @test !isimplicit(TableauRK4())
    @test !isdiagnonallyimplicit(TableauRK4())
    @test !isfullyimplicit(TableauRK4())
    
    @test typeof(TableauRK438()) <: Tableau
    @test order(TableauRK438()) == 4
    @test nstages(TableauRK438()) == 4

    @test typeof(TableauSSPRK2()) <: Tableau
    @test order(TableauSSPRK2()) == 2
    @test nstages(TableauSSPRK2()) == 2

    @test typeof(TableauSSPRK3()) <: Tableau
    @test order(TableauSSPRK3()) == 3
    @test nstages(TableauSSPRK3()) == 3

end


@testset "$(rpad("DIRK Tableaus",80))" begin

    @test typeof(TableauCrankNicolson()) <: Tableau
    @test order(TableauCrankNicolson()) == 2
    @test nstages(TableauCrankNicolson()) == 2
    @test issymmetric(TableauCrankNicolson())

    @test !isexplicit(TableauCrankNicolson())
    @test  isimplicit(TableauCrankNicolson())
    @test  isdiagnonallyimplicit(TableauCrankNicolson())
    @test !isfullyimplicit(TableauCrankNicolson())
    
    @test typeof(TableauCrouzeix()) <: Tableau
    @test order(TableauCrouzeix()) == 3
    @test nstages(TableauCrouzeix()) == 2

    @test !isexplicit(TableauCrouzeix())
    @test  isimplicit(TableauCrouzeix())
    @test  isdiagnonallyimplicit(TableauCrouzeix())
    @test !isfullyimplicit(TableauCrouzeix())
    
    @test typeof(TableauKraaijevangerSpijker()) <: Tableau
    @test order(TableauKraaijevangerSpijker()) == 2
    @test nstages(TableauKraaijevangerSpijker()) == 2

    @test typeof(TableauQinZhang()) <: Tableau
    @test order(TableauQinZhang()) == 2
    @test nstages(TableauQinZhang()) == 2
    @test issymplectic(TableauQinZhang())

end


@testset "$(rpad("FIRK Tableaus",80))" begin

    @test typeof(TableauImplicitEuler()) <: Tableau
    @test order(TableauImplicitEuler()) == 1
    @test nstages(TableauImplicitEuler()) == 1

    @test !isexplicit(TableauImplicitEuler())
    @test  isimplicit(TableauImplicitEuler())
    @test !isdiagnonallyimplicit(TableauImplicitEuler())
    @test  isfullyimplicit(TableauImplicitEuler())

    @test TableauBackwardEuler() == TableauImplicitEuler()

    @test typeof(TableauImplicitMidpoint()) <: Tableau
    @test order(TableauImplicitMidpoint()) == 2
    @test nstages(TableauImplicitMidpoint()) == 1
    @test issymplectic(TableauImplicitMidpoint())
    @test issymmetric(TableauImplicitMidpoint())

    @test !isexplicit(TableauImplicitMidpoint())
    @test  isimplicit(TableauImplicitMidpoint())
    @test !isdiagnonallyimplicit(TableauImplicitMidpoint())
    @test  isfullyimplicit(TableauImplicitMidpoint())

    @test typeof(TableauSRK3()) <: Tableau
    @test order(TableauSRK3()) == 4
    @test nstages(TableauSRK3()) == 3
    @test issymplectic(TableauSRK3())
    @test issymmetric(TableauSRK3())
    
end


@testset "$(rpad("PRK Tableaus",80))" begin

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
