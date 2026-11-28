# Rooted trees for the order conditions of a Runge-Kutta tableau (Butcher's theory).
#
# The trees of all orders up to the largest one requested so far are kept in one enumeration, by
# increasing order. A tree with the subtrees t₁, …, tₘ stores their positions in that enumeration,
# in non-increasing order, so that each tree is enumerated exactly once.

"""
    RootedTree

A rooted tree, with its level sequence (the depths of its nodes in preorder, the root at depth 1),
the positions of its subtrees in the enumeration of all trees by increasing order, its density ``γ(t)`` and its symmetry ``σ(t)``.
"""
struct RootedTree
    levels::Vector{Int}
    children::Vector{Int}
    γ::Int
    σ::Int
end

Base.copy(t::RootedTree) = RootedTree(copy(t.levels), copy(t.children), t.γ, t.σ)

"The order ``|t|`` of the tree `t`, its number of nodes."
tree_order(t::RootedTree) = length(t.levels)

"The level sequence of the tree `t`: the depths of its nodes in preorder, the root at depth 1."
levels(t::RootedTree) = copy(t.levels)

"""
The density ``γ(t)`` of the tree `t`: ``|t|`` times the densities of its subtrees, so that the
order condition of `t` reads ``b^T ψ(t) = 1/γ(t)``.
"""
density(t::RootedTree) = t.γ

"""
The symmetry ``σ(t)`` of the tree `t`: the order of its automorphism group, the product of the
symmetries of its subtrees and of ``m!`` for each subtree that occurs ``m`` times.
"""
symmetry(t::RootedTree) = t.σ

# The trees of order 1, 2, … in one enumeration; _TREE_STARTS[n] is the position of the first tree
# of order n, and _TREE_STARTS[end] the position after the last tree.
const _TREES = RootedTree[]
const _TREE_STARTS = [1]
const _TREE_LOCK = ReentrantLock()

# The tree of order n with the subtrees at the positions `children`, in non-increasing order.
function _make_tree(n, children)
    levels = [1]
    γ = n
    σ = 1
    multiplicity = 0
    for (i, k) in enumerate(children)
        t = _TREES[k]
        append!(levels, t.levels .+ 1)
        γ *= t.γ
        σ *= t.σ
        multiplicity = (i > 1 && children[i - 1] == k) ? multiplicity + 1 : 1
        σ *= multiplicity
    end
    return RootedTree(levels, copy(children), γ, σ)
end

# Push each tree of order n whose subtrees take the positions `children` and then subtrees of
# total order `rest`, each at a position no larger than `largest`.
function _push_trees!(trees, n, rest, largest, children)
    if rest == 0
        push!(trees, _make_tree(n, children))
        return trees
    end
    for k in min(largest, _TREE_STARTS[rest + 1] - 1):-1:1
        push!(children, k)
        _push_trees!(trees, n, rest - tree_order(_TREES[k]), k, children)
        pop!(children)
    end
    return trees
end

# The trees of order 1 … p, by increasing order, in a new vector. The enumeration grows under a
# lock, and the trees are shared with the cache; a caller that hands a tree on copies it.
function _trees_to(p::Integer)
    lock(_TREE_LOCK) do
        while length(_TREE_STARTS) <= p
            n = length(_TREE_STARTS)
            append!(_TREES, _push_trees!(RootedTree[], n, n - 1, length(_TREES), Int[]))
            push!(_TREE_STARTS, length(_TREES) + 1)
        end
        return _TREES[1:(_TREE_STARTS[p + 1] - 1)]
    end
end

"""
    rooted_trees(n)

The rooted trees of order `n`, as new copies of the enumeration.
"""
rooted_trees(n::Integer) = [copy(t) for t in _trees_to(n) if tree_order(t) == n]

# The stage vectors ψ(t) of the elementary weights bᵀψ(t) of `trees`, which hold every subtree
# before the tree, and the same vectors |ψ|(t) of |A|. With the subtrees t₁, …, tₘ of t,
# ψ(t) = (A ψ(t₁)) ⋯ (A ψ(tₘ)), entrywise, and ψ of the tree with one node is 𝟙. So a leaf below
# the root contributes A𝟙, never the stored nodes c.
function _stage_vectors(tab::Tableau{T}, trees) where {T}
    A = Matrix{T}(tab.a)
    absA = abs.(A)
    ψ = Vector{Vector{T}}(undef, length(trees))
    aψ = Vector{Vector{T}}(undef, length(trees))
    Aψ = Vector{Vector{T}}(undef, length(trees))
    absAψ = Vector{Vector{T}}(undef, length(trees))
    for (i, t) in enumerate(trees)
        v = ones(T, tab.s)
        w = ones(T, tab.s)
        for k in t.children
            v = v .* Aψ[k]
            w = w .* absAψ[k]
        end
        ψ[i], aψ[i] = v, w
        Aψ[i], absAψ[i] = A * v, absA * w
    end
    return ψ, aψ
end

"""
    stage_vectors(tab, p)

The stage vectors ``ψ(t)`` and ``|ψ|(t)`` of `tab` for the trees of order 1 … `p`, in the order of
their enumeration: the elementary weight of a tree ``t`` is ``b^T ψ(t)``, and ``|ψ|(t)`` is
``ψ(t)`` of ``|A|``.
"""
stage_vectors(tab::Tableau, p::Integer) = _stage_vectors(tab, _trees_to(p))
