using RungeKutta

@testset "$(rpad("Tableau List",80))" begin

    tab_list = TableauList

    for tab in tab_list
        @test typeof(tab()) <: Tableau
    end

end
