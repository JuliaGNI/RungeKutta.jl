using RungeKutta
using Test

import RootedTrees

const QQ = Rational{BigInt}

@testset "$(rpad("Rooted trees: counts, symmetry and density",80))" begin
    # the number of rooted trees of order n (OEIS A000081)
    counts = [1, 1, 2, 4, 9, 20, 48, 115, 286, 719]

    for n in 1:10
        trees = RungeKutta.rooted_trees(n)
        @test length(trees) == counts[n]
        @test all(t -> RungeKutta.tree_order(t) == n, trees)

        # every tree once: the canonical level sequences are distinct
        @test allunique(
            RootedTrees.canonical_representation(RootedTrees.rootedtree(RungeKutta.levels(t))).level_sequence
        for t in trees)

        # Σ 1/σ(t) = n^(n-1)/n!, the number of labelled rooted trees over n!
        @test sum(t -> 1 // QQ(RungeKutta.symmetry(t)), trees) ==
              QQ(n)^(n - 1) // factorial(big(n))

        # Σ 1/(σ(t) γ(t)) = 1/n, the coefficient of hⁿ in the Taylor series of exp(h)
        @test sum(t -> 1 // QQ(RungeKutta.symmetry(t) * RungeKutta.density(t)), trees) ==
              1 // QQ(n)

        @test RungeKutta.density(first(trees)) isa Int
        @test RungeKutta.symmetry(first(trees)) isa Int
    end

    @test isempty(RungeKutta.rooted_trees(0))

    # the level sequence starts at the root, at level 1, and lists the depths in preorder
    @test RungeKutta.levels(only(RungeKutta.rooted_trees(1))) == [1]
    @test sort([RungeKutta.levels(t) for t in RungeKutta.rooted_trees(3)]) ==
          [[1, 2, 2], [1, 2, 3]]
end

@testset "$(rpad("Rooted trees: no caller gets the cached storage",80))" begin
    # the level sequence is a copy
    t = first(RungeKutta.rooted_trees(4))
    l = RungeKutta.levels(t)
    l[end] = 100
    @test RungeKutta.levels(t) != l
    @test maximum(RungeKutta.levels(t)) <= 4

    # the trees are copies of the enumeration
    t.levels[end] = 100
    t.children[end] = 100
    @test maximum(RungeKutta.levels(first(RungeKutta.rooted_trees(4)))) <= 4
    @test maximum(first(RungeKutta.rooted_trees(4)).children) < 100
end

@testset "$(rpad("Rooted trees: RootedTrees.jl as an oracle",80))" begin
    tableaus = (
        TableauRK4(QQ),
        TableauHeun2(QQ),
        Tableau{QQ}(:IIIB2, 2, [1//2 0; 1//2 0], [1//2, 1//2], [0//1, 1//1]),
        # Lobatto IIIA(3), with a negative coefficient, so that |ψ| differs from ψ
        Tableau{QQ}(:IIIA3, 4, [0 0 0; 5//24 1//3 -1//24; 1//6 2//3 1//6],
            [1 // 6, 2 // 3, 1 // 6], [0, 1 // 2, 1]),
        TableauGauss(BigFloat, 3)
    )

    pmax = 8
    trees = reduce(vcat, [RungeKutta.rooted_trees(n) for n in 1:pmax])

    # the per-order counts and the trees themselves agree with RootedTrees.jl
    for n in 1:pmax
        ours = Set(RootedTrees.canonical_representation(RootedTrees.rootedtree(RungeKutta.levels(t))).level_sequence
        for t in trees if RungeKutta.tree_order(t) == n)
        theirs = Set(copy(t.level_sequence) for t in RootedTrees.RootedTreeIterator(n))
        @test length(ours) == count(t -> RungeKutta.tree_order(t) == n, trees)
        @test length(theirs) == length(ours)
        @test ours == theirs
    end

    for t in trees
        rt = RootedTrees.rootedtree(RungeKutta.levels(t))
        @test RungeKutta.density(t) == RootedTrees.density(rt)
        @test RungeKutta.symmetry(t) == RootedTrees.symmetry(rt)
    end

    for tab in tableaus
        T = eltype(tab)
        A = Matrix(tab.a)
        b = Vector(tab.b)
        A1 = A * ones(T, tab.s)
        ψ, aψ = RungeKutta.stage_vectors(tab, pmax)
        @test length(ψ) == length(aψ) == length(trees)
        @test eltype(eltype(ψ)) == eltype(eltype(aψ)) == T

        for (i, t) in enumerate(trees)
            Φ = sum(b .* ψ[i])
            Φref = RootedTrees.elementary_weight(RootedTrees.rootedtree(RungeKutta.levels(t)), A, b, A1)
            if T == QQ
                @test Φ isa QQ
                @test Φ == Φref
                # |ψ| is ψ of |A|
                @test sum(abs.(b) .* aψ[i]) ==
                      RootedTrees.elementary_weight(
                    RootedTrees.rootedtree(RungeKutta.levels(t)),
                    abs.(A), abs.(b), abs.(A) * ones(T, tab.s))
            else
                @test Φ isa BigFloat
                bound = 16 * RungeKutta.tree_order(t) * eps(BigFloat) *
                        sum(abs.(b) .* aψ[i])
                @test abs(Φ - Φref) <= bound
            end
        end
    end
end
