@testset "$(rpad("Tableau",80))" begin

    using RungeKutta: name, order, nstages, coefficients, weights, nodes

    for s in 1:5
        for T ∈ (Float64, BigFloat)

            a = rand(T,s,s)
            b = rand(T,s)
            c = rand(T,s)

            tab1 = Tableau(:Test, 2s, a, b, c)
            tab2 = Tableau(:Test, 2s, s, a, b, c)
            tab3 = Tableau{T}(:Test, 2s, a, b, c)
            tab4 = Tableau{T}(:Test, 2s, s, a, b, c)
            tabϵ = Tableau(:Test, 2s, a .+ 4eps(T) .* rand(T,s,s),
                                      b .+ 4eps(T) .* rand(T,s),
                                      c .+ 4eps(T) .* rand(T,s))

            @test tab1 ≈ tab2 ≈ tab3 ≈ tab4
            @test tab1 == tab2 == tab3 == tab4
            @test tab1 === tab2 === tab3 === tab4
            @test hash(tab1) == hash(tab2) == hash(tab3) == hash(tab4)

            @test tab1 ≈ tabϵ
            @test tab1 != tabϵ
            @test tab1 !== tabϵ
            @test hash(tab1) != hash(tabϵ)

            @test eltype(tab1) == eltype(tab2) == eltype(tab3) == eltype(tab4) == T

            @test name(tab1) == :Test
            @test order(tab1) == 2s
            @test nstages(tab1) == s
            @test coefficients(tab1) == a
            @test weights(tab1) == b
            @test nodes(tab1) == c

        end
    end

end
