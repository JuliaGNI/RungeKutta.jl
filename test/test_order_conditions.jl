@testset "$(rpad("Order Conditions",80))" begin

    for s in 1:5
        g = TableauGauss(s)

        for σ in 1:s
            @test check_order_conditions_b(g,σ) == true
            @test check_order_conditions_c(g,σ) == Array{Bool}(ones(s))
            @test check_order_conditions_d(g,σ) == Array{Bool}(ones(s))
        end

        @test satisfies_simplifying_assumption_b(g)
        @test satisfies_simplifying_assumption_c(g)
        @test satisfies_simplifying_assumption_d(g)
    end

    for s in 2:5
        @test  satisfies_simplifying_assumption_b(TableauLobattoIIIA(s))
        @test  satisfies_simplifying_assumption_b(TableauLobattoIIIB(s))
        @test  satisfies_simplifying_assumption_b(TableauLobattoIIIC(s))
        @test  satisfies_simplifying_assumption_b(TableauLobattoIIIC̄(s))
        @test  satisfies_simplifying_assumption_b(TableauLobattoIIID(s))
        @test  satisfies_simplifying_assumption_b(TableauLobattoIIIE(s))
        @test  satisfies_simplifying_assumption_b(TableauRadauIA(s))
        @test  satisfies_simplifying_assumption_b(TableauRadauIIA(s))

        @test  satisfies_simplifying_assumption_c(TableauLobattoIIIA(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIB(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIC(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIC̄(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIID(s))
        @test !satisfies_simplifying_assumption_c(TableauLobattoIIIE(s))
        @test !satisfies_simplifying_assumption_c(TableauRadauIA(s))
        @test  satisfies_simplifying_assumption_c(TableauRadauIIA(s))

        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIA(s))
        @test  satisfies_simplifying_assumption_d(TableauLobattoIIIB(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIC(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIC̄(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIID(s))
        @test !satisfies_simplifying_assumption_d(TableauLobattoIIIE(s))
        @test  satisfies_simplifying_assumption_d(TableauRadauIA(s))
        @test !satisfies_simplifying_assumption_d(TableauRadauIIA(s))
    end

end
