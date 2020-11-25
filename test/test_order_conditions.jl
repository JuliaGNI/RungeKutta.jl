@testset "$(rpad("Order Conditions",80))" begin

    s = 2
    g = TableauGauss(s)

    @test check_order_conditions_b(g,1) == true
    @test check_order_conditions_b(g,2) == true
    @test check_order_conditions_c(g,1) == Array{Bool}(ones(s))
    @test check_order_conditions_c(g,2) == Array{Bool}(ones(s))
    @test check_order_conditions_d(g,1) == Array{Bool}(ones(s))
    @test check_order_conditions_d(g,2) == Array{Bool}(ones(s))

end
