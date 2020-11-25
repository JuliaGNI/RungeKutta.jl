@testset "$(rpad("Tableau",80))" begin

    s = 2
    g = TableauGauss(s)

    @test check_symmetry(g) == Array{Bool}(ones(s,s))

end
