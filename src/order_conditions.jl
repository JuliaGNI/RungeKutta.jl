
# The default tolerance of the order-condition, symmetry and symplecticity checks: exactly zero
# for an exact element type, so that the check compares exactly, and 16 machine epsilons of the
# element type otherwise.
_default_tolerance(::Type{T}) where {T} = 16eps(T)
_default_tolerance(::Type{T}) where {T <: Union{Integer, Rational}} = zero(T)

# The right-hand side x / k of an order condition. For an integer or rational x it is the exact
# rational x // k, because x / k converts an integer x to Float64, where integers above 2^53
# that differ are equal. It is formed in BigInt, so that isapprox subtracts it from the left-hand
# side without the OverflowError of checked Int64 rational arithmetic.
_divide(x, k) = x / k
_divide(x::Union{Integer, Rational}, k::Integer) = big(x) // k

# The two sides of the simplifying assumptions B(k), C(k) and D(k) of a tableau: one pair for B,
# one pair per stage i for C and D. The residuals are the differences of the sides, with no
# tolerance; the checks compare the sides with `isapprox`.
function _sides_b(tab::Tableau{T}, k) where {T}
    b, c, s = tab.b, tab.c, tab.s
    return (mapreduce(i -> b[i] * c[i]^(k-1), +, 1:s), _divide(one(T), k))
end

function _sides_c(tab::Tableau, k)
    a, c, s = tab.a, tab.c, tab.s
    return [(mapreduce(j -> a[i, j] * c[j]^(k-1), +, 1:s), _divide(c[i]^k, k)) for i in 1:s]
end

function _sides_d(tab::Tableau, k)
    a, b, c, s = tab.a, tab.b, tab.c, tab.s
    return [(mapreduce(j -> b[j] * c[j]^(k-1) * a[j, i], +, 1:s),
                _divide(b[i] * (1 - c[i]^k), k)) for i in 1:s]
end

"The residual ``\\sum_j b_j c_j^{k-1} - 1/k`` of the simplifying assumption ``B(k)``."
residual_b(tab::Tableau, k) = -(_sides_b(tab, k)...)

"The residuals ``\\sum_j a_{ij} c_j^{k-1} - c_i^k/k``, ``i = 1, …, s``, of ``C(k)``."
residuals_c(tab::Tableau, k) = [l - r for (l, r) in _sides_c(tab, k)]

"The residuals ``\\sum_j b_j c_j^{k-1} a_{ji} - b_i (1 - c_i^k)/k``, ``i = 1, …, s``, of ``D(k)``."
residuals_d(tab::Tableau, k) = [l - r for (l, r) in _sides_d(tab, k)]

function check_order_conditions_b(tab::Tableau{T}, k; atol = _default_tolerance(T),
        rtol = _default_tolerance(T)) where {T}
    isapprox(_sides_b(tab, k)...; atol = atol, rtol = rtol)
end

function satisfies_simplifying_assumption_b(tab::Tableau, σ = tab.s; kwargs...)
    all([check_order_conditions_b(tab, k; kwargs...) for k in 1:σ])
end

function check_order_conditions_c(tab::Tableau{T}, k; atol = _default_tolerance(T),
        rtol = _default_tolerance(T)) where {T}
    [isapprox(l, r; atol = atol, rtol = rtol) for (l, r) in _sides_c(tab, k)]
end

function satisfies_simplifying_assumption_c(tab::Tableau, σ = tab.s; kwargs...)
    all(hcat([check_order_conditions_c(tab, k; kwargs...) for k in 1:σ]...))
end

function check_order_conditions_d(tab::Tableau{T}, k; atol = _default_tolerance(T),
        rtol = _default_tolerance(T)) where {T}
    [isapprox(l, r; atol = atol, rtol = rtol) for (l, r) in _sides_d(tab, k)]
end

function satisfies_simplifying_assumption_d(tab::Tableau, σ = tab.s; kwargs...)
    all(hcat([check_order_conditions_d(tab, k; kwargs...) for k in 1:σ]...))
end

@doc raw"""
Compute the weights by solving the simplifying assumption $B(s)$:
```math
\sum \limits_{j=1}^{s} b_{j} c_{j}^{k-1} = \frac{1}{k}  \qquad k = 1 , \, ... , \, s .
```
"""
function solve_simplifying_assumption_b(c::AbstractVector{T}) where {T}
    s = length(c)
    M = [c[j]^(k-1) for k in 1:s, j in 1:s]
    r = [1 / T(k) for k in 1:s]
    M \ r
end

@doc raw"""
Compute the coefficients by solving the simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function solve_simplifying_assumption_c(c::AbstractVector{T}) where {T}
    s = length(c)
    M = [c[j]^(k-1) for k in 1:s, j in 1:s]

    row(i) = begin
        r = [c[i]^k / k for k in 1:s]
        M \ r
    end

    vcat([transpose(row(i)) for i in 1:s]...)
end

@doc raw"""
Compute the coefficients by solving the simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function solve_simplifying_assumption_d(b::AbstractVector{T}, c::AbstractVector{T}) where {T}
    @assert axes(b) == axes(c)

    s = length(c)
    M = [b[i] * c[i]^(k-1) for k in 1:s, i in 1:s]

    row(j) = begin
        r = [b[j] / k * (1 - c[j]^k) for k in 1:s]
        M \ r
    end

    hcat([row(j) for j in 1:s]...)
end

function solve_simplifying_assumption_d(b::AbstractVector, c::AbstractVector)
    solve_simplifying_assumption_d(promote(b, c)...)
end

# The residual bᵀψ(t) - 1/γ(t) of the order condition of the tree t.
_tree_residual(b, ψ, t::RootedTree) = sum(b .* ψ) - _divide(one(eltype(b)), density(t))

# The number of leading k = 1, …, n for which f(k) holds.
_leading_count(f, n) = something(findfirst(!f, 1:n), n + 1) - 1

# The classical order that the simplifying assumptions B(q_B), C(η) and D(ζ) give by Butcher's
# theorem: B(p), C(η) and D(ζ) with p ≤ η + ζ + 1 and p ≤ 2η + 2 give order p.
_order_bound(qB, η, ζ) = min(qB, η + ζ + 1, 2η + 2)

@doc raw"""
    order_conditions(tab::Tableau{T}, p)

Return the order conditions of `tab` for the rooted trees of order 1 to `p`: a vector with one
entry `(tree = t, residual = r)` per tree ``t``, by increasing order of the trees, with

```math
r = b^T ψ(t) - \frac{1}{γ(t)} ,
```

the elementary weight ``b^T ψ(t)`` minus the inverse of the density ``γ(t)`` of ``t``. The stage
vector ``ψ(t)`` is built from ``A`` alone: the tree with one node has ``ψ = 𝟙``, and a tree with
the subtrees ``t_1, …, t_m`` has ``ψ(t) = (A ψ(t_1)) ⋯ (A ψ(t_m))``, entrywise. So the conditions
use the row sums ``A 𝟙``, not the stored nodes ``c``; [`isconsistent`](@ref) tells whether the
two agree. The classical order of the tableau is at least ``p`` where every residual of the trees
of order ``p`` or less is zero.

No tolerance is applied, and no residual is simplified. A residual has the element type `T` of
the tableau, such as `Float64`, `BigFloat`, `Rational{BigInt}` or a symbolic type, and
`Rational{BigInt}` for an integer `T`. `p = 0` gives an empty vector.

Throws an `ArgumentError` when `p` is negative or when a coefficient of `a` or `b` is not finite.
"""
function order_conditions(tab::Tableau{T}, p::Integer) where {T}
    p >= 0 || throw(ArgumentError("the order p must be nonnegative"))
    _check_finite(tab)
    b = Vector{T}(tab.b)
    trees = _trees_to(p)
    ψ, _ = _stage_vectors(tab, trees)
    return [(tree = copy(t), residual = _tree_residual(b, ψ[i], t))
            for (i, t) in enumerate(trees)]
end

@doc raw"""
    computed_order(tab::Tableau{T}; pmax = 10, rtol)

Return the classical order of `tab`, computed from its order conditions, as a named tuple
`(order = p, route, exact)`, with `p` an `Int`.

1. For each order ``n`` from 1 to `pmax`, the order conditions of the trees of order ``n`` are
   tested (see [`order_conditions`](@ref)). Where one fails, the result is
   `(order = n - 1, route = :trees, exact = true)`.
2. Where all hold, the simplifying assumptions ``B(q_B)``, ``C(η)`` and ``D(ζ)`` are counted on
   the nodes ``c = A 𝟙`` (not on the stored `c`), with ``q_B ≤ 2s`` and ``η, ζ ≤ s``, and give
   the order ``p_{BCD} = \min(q_B, η + ζ + 1, 2η + 2)`` by Butcher's theorem.
3. The result is `order = max(pmax, p_BCD)`, with `route = :simplifying_assumptions` where
   ``p_{BCD} >`` `pmax` and `route = :trees` otherwise, and `exact = order ≥ q_B`. Where `exact`
   is `false`, the order is a lower bound: a larger `pmax` can give a larger order.

The order condition of a tree ``t`` holds when

```math
|b^T ψ(t) - 1/γ(t)| ≤ \mathrm{rtol} \cdot |t| \cdot |b|^T |ψ|(t) ,
```

where ``|ψ|(t)`` is ``ψ(t)`` built from ``|A|``: the first-order bound on the change of the
elementary weight when every entry of ``A`` and ``b`` changes by a relative `rtol`. For the tall
trees it is the bound of [`linear_order`](@ref). The simplifying assumptions are checked with
[`check_order_conditions_b`](@ref), [`check_order_conditions_c`](@ref) and
[`check_order_conditions_d`](@ref), with `atol = rtol = rtol`. The default `rtol` is
`16 * eps(T)` for a floating-point `T` and `0` for a rational one, where every test is exact.

The classical order is at most the [`linear_order`](@ref), the order of the stability function,
and can be lower. It does not read the stated order `order(tab)`.

Throws an `ArgumentError` when `pmax < 1` or when a coefficient of `a` or `b` is not finite.
"""
function computed_order(tab::Tableau{T}; pmax::Integer = 10,
        rtol = _default_rtol(T)) where {T <: Real}
    pmax >= 1 || throw(ArgumentError("pmax must be at least 1"))
    _check_finite(tab)
    b = Vector{T}(tab.b)
    absb = abs.(b)
    trees = _trees_to(pmax)
    ψ, aψ = _stage_vectors(tab, trees)
    for (i, t) in enumerate(trees)
        abs(_tree_residual(b, ψ[i], t)) <= rtol * tree_order(t) * sum(absb .* aψ[i]) ||
            return (order = tree_order(t) - 1, route = :trees, exact = true)
    end

    s = tab.s
    tabA1 = Tableau{T}(tab.name, 0, s, tab.a, tab.b, tab.a * ones(T, s))
    qB = _leading_count(k -> check_order_conditions_b(tabA1, k; atol = rtol, rtol = rtol), 2s)
    η = _leading_count(k -> all(check_order_conditions_c(tabA1, k; atol = rtol, rtol = rtol)), s)
    ζ = _leading_count(k -> all(check_order_conditions_d(tabA1, k; atol = rtol, rtol = rtol)), s)
    pBCD = _order_bound(qB, η, ζ)
    p = max(pmax, pBCD)
    return (
        order = p, route = pBCD > pmax ? :simplifying_assumptions : :trees, exact = p >= qB)
end

@doc raw"""
    stage_order(tab::Tableau; kwargs...)

Return the stage order of `tab`, as an `Int`: the largest ``q ≤ s`` such that the simplifying
assumptions ``C(1), …, C(q)`` hold,

```math
\sum_{j=1}^{s} a_{ij} c_j^{k-1} = \frac{c_i^k}{k} , \qquad i = 1, …, s , \; k = 1, …, q ,
```

on the stored nodes ``c``, as [`check_order_conditions_c`](@ref) tests them, with its defaults and
the keywords `kwargs` (`atol`, `rtol`). It is 0 where ``c ≠ A 𝟙`` (see [`isconsistent`](@ref)).

Throws an `ArgumentError` when a coefficient of `a` or `b` is not finite.
"""
function stage_order(tab::Tableau; kwargs...)
    _check_finite(tab)
    return _leading_count(k -> all(check_order_conditions_c(tab, k; kwargs...)), tab.s)
end

@doc raw"""
    isconsistent(tab::Tableau; kwargs...)

Return `true` where the stored nodes of `tab` are the row sums of its coefficients, ``c = A 𝟙``,
that is where ``C(1)`` holds as [`check_order_conditions_c`](@ref) tests it, with its defaults and
the keywords `kwargs` (`atol`, `rtol`). The order conditions of [`order_conditions`](@ref) and
[`computed_order`](@ref) use ``A 𝟙`` in either case.

Throws an `ArgumentError` when a coefficient of `a` or `b` is not finite.
"""
function isconsistent(tab::Tableau; kwargs...)
    _check_finite(tab)
    return all(check_order_conditions_c(tab, 1; kwargs...))
end

@doc raw"""
    error_constant(tab::Tableau{T}; pmax = 10, rtol)

Return the error constant of `tab`, of type `float(T)`: the 2-norm of

```math
\frac{1}{σ(t)} \left( b^T ψ(t) - \frac{1}{γ(t)} \right)
```

over the rooted trees ``t`` of order ``p + 1``, with ``σ(t)`` the symmetry of ``t`` and ``p`` the
classical order from [`computed_order`](@ref) with `pmax` and `rtol`. These are the coefficients
``e(t)`` of the leading term ``h^{p+1} \sum_{|t| = p+1} e(t) F(t)(y_0)`` of the local error, with
``F(t)`` the elementary differential of ``t``.

Throws an `ArgumentError` when the computed order is a lower bound (`exact` is `false`; a larger
`pmax` can help), when `pmax < 1`, or when a coefficient of `a` or `b` is not finite.
"""
function error_constant(tab::Tableau{T}; pmax::Integer = 10,
        rtol = _default_rtol(T)) where {T <: Real}
    _check_finite(tab)
    r = computed_order(tab; pmax = pmax, rtol = rtol)
    r.exact || throw(ArgumentError("the computed order is a lower bound; increase pmax"))
    b = Vector{T}(tab.b)
    trees = _trees_to(r.order + 1)
    ψ, _ = _stage_vectors(tab, trees)
    F = float(T)
    e = [F(_tree_residual(b, ψ[i], t) / symmetry(t))
         for (i, t) in enumerate(trees) if tree_order(t) == r.order + 1]
    return LinearAlgebra.norm(e)
end
