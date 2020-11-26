using RungeKutta: name, order, nstages, coefficients, weights, nodes

@testset "$(rpad("ERK Tableaus",80))" begin

    @test typeof(TableauExplicitEuler()) <: Tableau
    @test order(TableauExplicitEuler()) == 1
    @test nstages(TableauExplicitEuler()) == 1

    @test TableauExplicitEuler() == TableauExplicitEuler()

    @test typeof(TableauExplicitMidpoint()) <: Tableau
    @test order(TableauExplicitMidpoint()) == 2
    @test nstages(TableauExplicitMidpoint()) == 2

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

    @test typeof(TableauRK438()) <: Tableau
    @test order(TableauRK438()) == 4
    @test nstages(TableauRK438()) == 4

    @test typeof(TableauSSPRK3()) <: Tableau
    @test order(TableauSSPRK3()) == 3
    @test nstages(TableauSSPRK3()) == 3

end


@testset "$(rpad("DIRK Tableaus",80))" begin

    @test typeof(TableauKraaijevangerSpijker()) <: Tableau
    @test order(TableauKraaijevangerSpijker()) == 2
    @test nstages(TableauKraaijevangerSpijker()) == 2

    @test typeof(TableauQinZhang()) <: Tableau
    @test order(TableauQinZhang()) == 2
    @test nstages(TableauQinZhang()) == 2
    @test issymplectic(TableauQinZhang())

    @test typeof(TableauCrouzeix()) <: Tableau
    @test order(TableauCrouzeix()) == 3
    @test nstages(TableauCrouzeix()) == 2

end


@testset "$(rpad("FIRK Tableaus",80))" begin

    @test typeof(TableauImplicitEuler()) <: Tableau
    @test order(TableauImplicitEuler()) == 1
    @test nstages(TableauImplicitEuler()) == 1

    @test TableauBackwardEuler() == TableauImplicitEuler()

    @test typeof(TableauImplicitMidpoint()) <: Tableau
    @test order(TableauImplicitMidpoint()) == 2
    @test nstages(TableauImplicitMidpoint()) == 1
    @test issymplectic(TableauImplicitMidpoint())
    @test issymmetric(TableauImplicitMidpoint())

    @test typeof(TableauSRK3()) <: Tableau
    @test order(TableauSRK3()) == 4
    @test nstages(TableauSRK3()) == 3
    @test issymplectic(TableauSRK3())
    @test issymmetric(TableauSRK3())

    @test typeof(TableauCrankNicolson()) <: Tableau
    @test order(TableauCrankNicolson()) == 2
    @test nstages(TableauCrankNicolson()) == 2
    @test issymmetric(TableauCrankNicolson())
    
end
