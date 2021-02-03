using RungeKutta: name, order, eachstage, nstages

@testset "$(rpad("Partitioned Tableau",80))" begin

    for s in 1:5
        for T ∈ (Float64, BigFloat)

            q = Tableau(:qTest, 2s, rand(T,s,s), rand(T,s), rand(T,s))
            p = Tableau(:pTest, 2s, rand(T,s,s), rand(T,s), rand(T,s))
            
            tab1 = PartitionedTableau{T}(:Test, 2s, s, q, p)
            tab2 = PartitionedTableau{T}(:Test, q, p)
            tab3 = PartitionedTableau(:Test, q, p)
            tab4 = PartitionedTableau(:Test, q)

            @test tab1 ≈ tab2 ≈ tab3
            @test tab1 == tab2 == tab3 != tab4
            @test hash(tab1) == hash(tab2) == hash(tab3) != hash(tab4)

            @test isequal(tab1, tab2)
            @test isequal(tab1, tab3)
            @test !isequal(tab1, tab4)

            @test eltype(tab1) == eltype(tab2) == eltype(tab3) == eltype(tab4) == T

            @test name(tab1) == :Test
            @test order(tab1) == 2s
            @test eachstage(tab1) == 1:s
            @test nstages(tab1) == s

        end

    end

end
