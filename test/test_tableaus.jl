using RungeKutta: name, order, nstages, coefficients, weights, nodes

@testset "$(rpad("ERK Tableaus",80))" begin

    @test typeof(TableauExplicitEuler()) <: Tableau
    @test order(TableauExplicitEuler()) == 1
    @test nstages(TableauExplicitEuler()) == 1

    @test typeof(TableauExplicitMidpoint()) <: Tableau
    @test order(TableauExplicitMidpoint()) == 2
    @test nstages(TableauExplicitMidpoint()) == 2

    @test typeof(TableauHeun()) <: Tableau
    @test order(TableauHeun()) == 2
    @test nstages(TableauHeun()) == 2

    @test typeof(TableauRunge()) <: Tableau
    @test order(TableauRunge()) == 2
    @test nstages(TableauRunge()) == 2

    @test typeof(TableauKutta()) <: Tableau
    @test order(TableauKutta()) == 3
    @test nstages(TableauKutta()) == 3

    @test typeof(TableauERK416()) <: Tableau
    @test order(TableauERK416()) == 4
    @test nstages(TableauERK416()) == 4

    @test typeof(TableauERK438()) <: Tableau
    @test order(TableauERK438()) == 4
    @test nstages(TableauERK438()) == 4

end


@testset "$(rpad("DIRK Tableaus",80))" begin

    @test typeof(TableauCrouzeix()) <: Tableau
    @test order(TableauCrouzeix()) == 3
    @test nstages(TableauCrouzeix()) == 2

end


@testset "$(rpad("FIRK Tableaus",80))" begin

    @test typeof(TableauImplicitEuler()) <: Tableau
    @test order(TableauImplicitEuler()) == 1
    @test nstages(TableauImplicitEuler()) == 1

    @test typeof(TableauImplicitMidpoint()) <: Tableau
    @test order(TableauImplicitMidpoint()) == 2
    @test nstages(TableauImplicitMidpoint()) == 1

    @test typeof(TableauSRK3()) <: Tableau
    @test order(TableauSRK3()) == 4
    @test nstages(TableauSRK3()) == 3

end
