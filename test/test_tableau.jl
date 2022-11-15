using RungeKutta: name, order, eachstage, nstages, coefficients, weights, nodes, to_array, to_file, from_file

@testset "$(rpad("Tableau",80))" begin

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
            @test hash(tab1) == hash(tab2) == hash(tab3) == hash(tab4)

            @test tab1 ≈ tabϵ
            @test tab1 != tabϵ
            @test hash(tab1) != hash(tabϵ)

            @test isequal(tab1, tab2)
            @test isequal(tab1, tab3)
            @test isequal(tab1, tab4)
            @test !isequal(tab1, tabϵ)

            @test eltype(tab1) == eltype(tab2) == eltype(tab3) == eltype(tab4) == T

            @test name(tab1) == :Test
            @test order(tab1) == 2s
            @test eachstage(tab1) == 1:s
            @test nstages(tab1) == s
            @test coefficients(tab1) == a
            @test weights(tab1) == b
            @test nodes(tab1) == c

            @test ismissing(tab1.R∞)
            @test ismissing(tab2.R∞)
            @test ismissing(tab3.R∞)
            @test ismissing(tab4.R∞)

            @test tab1.a == tab2.a == tab3.a == tab4.a == a
            @test tab1.b == tab2.b == tab3.b == tab4.b == b
            @test tab1.c == tab2.c == tab3.c == tab4.c == c

            @test tab1.â == tab2.â == tab3.â == tab4.â == zero(a)
            @test tab1.b̂ == tab2.b̂ == tab3.b̂ == tab4.b̂ == zero(b)
            @test tab1.ĉ == tab2.ĉ == tab3.ĉ == tab4.ĉ == zero(c)
            
            @test tab1 == Tableau(tab1.name, tab1.o, to_array(tab1))
            @test tab1 == Tableau(tab1.name, tab1.o, convert(Matrix, tab1))
            @test tab1 == Tableau(tab1.name, tab1.o, convert(Array, tab1))
            @test tab1 == convert(Tableau, convert(Matrix, tab1); name=tab1.name, o=tab1.o)
            @test tab1 == convert(Tableau, convert(Array,  tab1); name=tab1.name, o=tab1.o)
            @test convert(Matrix, tab1) == convert(Matrix{T}, tab1)
            @test convert(Array,  tab1) == convert(Array{T},  tab1)

            @test startswith(repr(tab1), "\nRunge-Kutta Tableau")
            @test startswith(repr(MIME("text/markdown"), tab1), "Runge-Kutta Tableau")

        end


        T1 = BigFloat
        T2 = Float64

        a = rand(T1,s,s)
        b = rand(T1,s)
        c = rand(T1,s)

        tab1 = Tableau{T1}(:Test, 2s, s, a, b, c)
        tab2 = Tableau{T2}(:Test, 2s, s, a, b, c)

        @test tab2.â == T2.(tab1.a .- tab2.a)
        @test tab2.b̂ == T2.(tab1.b .- tab2.b)
        @test tab2.ĉ == T2.(tab1.c .- tab2.c)

        tmp = mktempdir()
        to_file(tmp, tab2)
        tabf = from_file(tmp, "Test")
        rm(tmp, recursive=true)

        @test tabf.o == tab2.o
        @test tabf.s == tab2.s
        @test tabf.a == tab2.a
        @test tabf.b == tab2.b
        @test tabf.c == tab2.c

    end

end
