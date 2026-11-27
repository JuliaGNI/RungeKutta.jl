# Polynomials are vectors of coefficients in ascending powers: `p[k + 1]` is the coefficient of
# `z^k`. The decisions on the real line (signs, roots) are made in `Rational{BigInt}`, into which
# every coefficient of a `Float32`, `Float64`, `BigFloat` or rational tableau converts exactly.

const _Exact = Rational{BigInt}

"""
    StabilityFunction{T}

The stability function ``R(z) = P(z) / Q(z)`` of a Runge-Kutta tableau, the factor by which one
step multiplies the solution of ``\\dot{y} = λ y`` with ``z = h λ``. The fields `P` and `Q` are
`Vector{T}`s of length ``s + 1``, the coefficients of ``P(z) = \\det(I - z A + z 𝟙 b^T)`` and
``Q(z) = \\det(I - z A)`` in ascending powers of ``z``; none of them is trimmed or rounded.

`R(z)` evaluates the quotient `P(z) / Q(z)` for a real or complex `z`, as computed, also at a pole.

See [`stability_function`](@ref).
"""
struct StabilityFunction{T}
    P::Vector{T}
    Q::Vector{T}
end

(R::StabilityFunction)(z) = evalpoly(z, R.P) / evalpoly(z, R.Q)

function _check_finite(tab::Tableau)
    all(isfinite, tab.a) && all(isfinite, tab.b) ||
        throw(ArgumentError("the tableau has a coefficient that is not finite"))
    return nothing
end

# The coefficients of det(I - z M) in ascending powers of z, by the Faddeev–LeVerrier recursion.
# It uses the ring operations and the division by an integer only, so it is exact in a rational
# or a symbolic element type.
function _det_coefficients(M::AbstractMatrix{T}) where {T}
    s = size(M, 1)
    c = zeros(T, s + 1)
    c[1] = one(T)
    N = zeros(T, s, s)
    for i in 1:s
        N[i, i] = one(T)
    end
    for k in 1:s
        N = M * N
        c[k + 1] = -sum(N[i, i] for i in 1:s) / k
        for i in 1:s
            N[i, i] += c[k + 1]
        end
    end
    return c
end

"""
    stability_function(tab::Tableau{T})

Return the [`StabilityFunction`](@ref) `R` of `tab`: the coefficient vectors `R.P` and `R.Q`, of
element type `T` and length ``s + 1``, of

```math
R(z) = 1 + z b^T (I - z A)^{-1} 𝟙 = \\frac{P(z)}{Q(z)} , \\qquad
P(z) = \\det(I - z A + z 𝟙 b^T) , \\qquad
Q(z) = \\det(I - z A) .
```

The coefficients come from the Faddeev–LeVerrier recursion, which uses the ring operations and
the division by an integer only. They are exact for a rational element type such as
`Rational{BigInt}` and for a symbolic one, and carry the round-off of that recursion for a
floating-point one. No coefficient is trimmed: a coefficient that is zero in exact arithmetic can
hold round-off, and the functions that read degrees count it as zero with their keyword `rtol`.

Throws an `ArgumentError` when a coefficient of `a` or `b` is not finite.
"""
function stability_function(tab::Tableau{T}) where {T}
    _check_finite(tab)
    A = Matrix{T}(tab.a)
    b = Vector{T}(tab.b)
    Q = _det_coefficients(A)
    P = _det_coefficients(A - ones(T, tab.s) * transpose(b))
    return StabilityFunction{T}(P, Q)
end

# The default relative tolerance of the zero count: none for an exact element type.
_default_rtol(::Type{T}) where {T <: AbstractFloat} = 16 * eps(T)
_default_rtol(::Type{<:Real}) = 0

# The coefficients of P and Q with each one of magnitude at most rtol · m set to zero, where m is
# the largest magnitude of a coefficient of P and Q; and m.
function _zero_counted(R::StabilityFunction, rtol)
    m = max(maximum(abs, R.P), maximum(abs, R.Q))
    return _count_zeros(R.P, rtol * m), _count_zeros(R.Q, rtol * m), m
end

_count_zeros(p, tol) = [abs(c) <= tol ? zero(c) : c for c in p]

# The degree of a polynomial whose constant coefficient is not zero.
_degree(p) = findlast(!iszero, p) - 1

"""
    stability_at_infinity(tab::Tableau{T}; rtol)

Return ``R(∞)``, the limit of the stability function ``R(z) = P(z)/Q(z)`` of `tab` as
``|z| → ∞``, of type `T`: the quotient of the leading coefficients of ``P`` and ``Q`` where their
degrees are equal, `zero(T)` where ``\\deg P < \\deg Q``, and the positive infinity of `T` where
``\\deg P > \\deg Q`` (`1//0` for a rational `T`). It is computed from the coefficients of
[`stability_function`](@ref), and does not read the field `R∞` of `tab`.

The degrees follow from a zero count: a coefficient ``c`` of ``P`` or ``Q`` counts as zero when
``|c| ≤ \\mathrm{rtol} \\cdot m``, with ``m`` the largest magnitude of a coefficient of ``P`` and
``Q``. The default `rtol` is `16 * eps(T)` for a floating-point `T` and `0` for a rational one.
"""
function stability_at_infinity(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    P, Q, _ = _zero_counted(stability_function(tab), rtol)
    dP, dQ = _degree(P), _degree(Q)
    dP < dQ && return zero(T)
    dP > dQ && return one(T) / zero(T)
    return P[dP + 1] / Q[dQ + 1]
end

"""
    linear_order(tab::Tableau{T}; rtol)

Return the linear order of `tab`, as an `Int`: the largest ``q`` with
``b^T A^{k-1} 𝟙 = 1/k!`` for ``k = 1, …, q``. It is the order of the stability function ``R(z)``
as an approximation of ``\\exp(z)``, and so the order of convergence on a linear problem with
constant coefficients; the classical order, from the rooted trees, can be lower. It is at most
``\\deg P + \\deg Q ≤ 2s``, because a rational ``R`` is not ``\\exp``, so ``k`` runs to ``2s``.

A residual counts as zero when ``|b^T A^{k-1} 𝟙 - 1/k!| ≤ \\mathrm{rtol} / k!``. The default
`rtol` is `16 * eps(T)` for a floating-point `T` and `0` for a rational one.
"""
function linear_order(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    _check_finite(tab)
    A = Matrix{T}(tab.a)
    b = Vector{T}(tab.b)
    v = ones(T, tab.s)
    kfactorial = one(T)
    q = 0
    for k in 1:(2 * tab.s)
        kfactorial *= k
        abs(sum(b .* v) - inv(kfactorial)) <= rtol * inv(kfactorial) || break
        q = k
        v = A * v
    end
    return q
end

"""
    pade_type(tab::Tableau{T}; rtol)

Return `(deg P, deg Q)` where the stability function ``R = P/Q`` of `tab` is the Padé
approximation of ``\\exp(z)`` of that type, that is where the [`linear_order`](@ref) of `tab`
equals ``\\deg P + \\deg Q``, and `nothing` otherwise. The degrees are `Int`s.

The degrees follow from the zero count of [`stability_at_infinity`](@ref), and the linear order
from its own residual test, both with the same `rtol`. The default `rtol` is `16 * eps(T)` for a
floating-point `T` and `0` for a rational one.
"""
function pade_type(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    P, Q, _ = _zero_counted(stability_function(tab), rtol)
    dP, dQ = _degree(P), _degree(Q)
    return linear_order(tab; rtol) == dP + dQ ? (dP, dQ) : nothing
end

# Polynomial arithmetic on ascending coefficient vectors. A zero polynomial may be empty.

function _padd(p, q)
    n = max(length(p), length(q))
    return [(k <= length(p) ? p[k] : zero(eltype(p))) +
            (k <= length(q) ? q[k] : zero(eltype(q)))
            for k in 1:n]
end

_psub(p, q) = _padd(p, -q)

function _pmul(p::AbstractVector{T}, q::AbstractVector{T}) where {T}
    (isempty(p) || isempty(q)) && return T[]
    r = zeros(T, length(p) + length(q) - 1)
    for i in eachindex(p), j in eachindex(q)

        r[i + j - 1] += p[i] * q[j]
    end
    return r
end

# p(-z)
_reflect(p) = [isodd(k) ? -p[k + 1] : p[k + 1] for k in 0:(length(p) - 1)]

# |p(iy)|² as a polynomial in y, for real coefficients: p(iy) = re(y) + i im(y).
function _abs2_on_imaginary_axis(p::AbstractVector{T}) where {T}
    re = [iseven(k) ? (iseven(k ÷ 2) ? p[k + 1] : -p[k + 1]) : zero(T)
          for k in 0:(length(p) - 1)]
    im = [isodd(k) ? (iseven(k ÷ 2) ? p[k + 1] : -p[k + 1]) : zero(T)
          for k in 0:(length(p) - 1)]
    return _padd(_pmul(re, re), _pmul(im, im))
end

_ptrim(p) = p[1:something(findlast(!iszero, p), 0)]

_pderivative(p) = [k * p[k + 1] for k in 1:(length(p) - 1)]

function _pdivrem(a::Vector{_Exact}, b::Vector{_Exact})
    a, b = _ptrim(a), _ptrim(b)
    q = zeros(_Exact, max(length(a) - length(b) + 1, 0))
    r = copy(a)
    for k in length(q):-1:1
        f = r[k + length(b) - 1] / b[end]
        q[k] = f
        for j in eachindex(b)
            r[k + j - 1] -= f * b[j]
        end
    end
    return q, _ptrim(r)
end

_pdiv(a, b) = first(_pdivrem(a, b))

function _pgcd(a::Vector{_Exact}, b::Vector{_Exact})
    a, b = _ptrim(a), _ptrim(b)
    while !isempty(b)
        a, b = b, last(_pdivrem(a, b))
    end
    return a ./ a[end]
end

# The product of the factors of odd multiplicity of p, by the square-free factorisation of Yun:
# p changes its sign exactly at the real roots of this product.
function _odd_multiplicity_part(p::Vector{_Exact})
    dp = _pderivative(p)
    g = _pgcd(p, dp)
    b = _pdiv(p, g)
    d = _psub(_pdiv(dp, g), _pderivative(b))
    odd = _Exact[1]
    multiplicity = 1
    while length(b) > 1
        a = _pgcd(b, d)
        isodd(multiplicity) && (odd = _pmul(odd, a))
        b = _pdiv(b, a)
        d = _psub(_pdiv(d, a), _pderivative(b))
        multiplicity += 1
    end
    return odd
end

# The Sturm sequence of a polynomial with simple roots, and its number of sign variations at x.
function _sturm_sequence(p::Vector{_Exact})
    seq = [p, _pderivative(p)]
    while length(seq[end]) > 1
        push!(seq, -last(_pdivrem(seq[end - 1], seq[end])))
    end
    return seq
end

function _sign_variations(seq, x)
    signs = filter(!iszero, [sign(evalpoly(x, p)) for p in seq if !isempty(p)])
    return count(k -> signs[k] != signs[k + 1], 1:(length(signs) - 1))
end

# The roots in (0, ∞) of p, which has simple roots and p(0) ≠ 0, in increasing order. Each root is
# bracketed by bisection on the Sturm count until the bracket holds no other root and both of its
# ends round to one value of F; the midpoint of the bracket is returned.
function _positive_roots(p::Vector{_Exact}, ::Type{F}) where {F}
    roots = _Exact[]
    length(p) == 1 && return roots
    seq = _sturm_sequence(p)
    bound = 1 + maximum(abs, p[1:(end - 1)]) / abs(p[end])
    vbound = _sign_variations(seq, bound)
    a = zero(_Exact)
    va = _sign_variations(seq, a)
    while va > vbound
        lo, hi, vhi = a, bound, vbound
        for _ in 1:(4 * precision(F) + 4096)
            F(lo) == F(hi) && vhi == va - 1 && break
            mid = (lo + hi) / 2
            vmid = _sign_variations(seq, mid)
            if vmid < va
                hi, vhi = mid, vmid
            else
                lo = mid
            end
        end
        push!(roots, (lo + hi) / 2)
        a, va = hi, vhi
    end
    return roots
end

# The coefficients of f, a polynomial that is zero at 0, in exact arithmetic, with each one of
# magnitude at most tol counted as zero and the factor xʲ of the lowest nonzero one divided out;
# `nothing` where f is negative just right of 0.
function _positive_part(f, tol)
    g = _Exact[abs(c) <= tol ? 0 : _Exact(c) for c in f]
    k = findfirst(!iszero, g)
    k === nothing && return _Exact[]
    g[k] < 0 && return nothing
    return _ptrim(g[k:end])
end

# The supremum of the r ≥ 0 with f ≥ 0 on [0, r], as F: the first root at which f changes its sign
# to negative, that is a root of odd multiplicity. Where `touch` is true, a sign change counts only
# where f, between that root and the next, falls below -tol Σₖ xᵏ, the bound of a change of each
# coefficient by tol: a point where f touches zero then does not end the interval when round-off
# splits it into two close roots.
function _nonnegative_until(f, tol, ::Type{F}; touch = false) where {F}
    g = _positive_part(f, tol)
    g === nothing && return zero(F)
    isempty(g) && return F(Inf)
    roots = _positive_roots(_odd_multiplicity_part(g), F)
    for (i, ρ) in enumerate(roots)
        if !touch
            return F(ρ)
        elseif i == length(roots)
            g[end] < 0 && return F(ρ)
        else
            x = (ρ + roots[i + 1]) / 2
            evalpoly(x, g) < -tol * evalpoly(x, ones(_Exact, length(g))) && return F(ρ)
        end
    end
    return F(Inf)
end

_isnonnegative(f, tol) = isinf(_nonnegative_until(f, tol, Float64))

# Whether every root of q, whose constant coefficient is not zero, has a positive real part: the
# Routh–Hurwitz test of q(-z), each entry of the first column of the Routh array of the sign of
# the leading coefficient.
function _has_roots_in_right_half_plane_only(q)
    c = reverse(_ptrim(_Exact.(_reflect(q))))
    r1, r2 = c[1:2:end], c[2:2:end]
    σ = sign(c[1])
    for _ in 1:(length(c) - 1)
        append!(r2, zeros(_Exact, length(r1) - length(r2)))
        sign(r2[1]) == σ || return false
        r1, r2 = r2, [r1[j + 1] - r1[1] * r2[j + 1] / r2[1] for j in 1:(length(r1) - 1)]
    end
    return true
end

"""
    isastable(tab::Tableau{T}; rtol) where {T <: Real}

Return whether `tab` is A-stable, that is whether ``|R(z)| ≤ 1`` for every ``z`` with
``\\mathrm{Re}\\, z ≤ 0``. It is decided from the coefficients of ``R = P/Q`` without sampling
``z``, by the two conditions that are together equivalent to it:

1. ``Q`` has no zero with ``\\mathrm{Re}\\, z ≤ 0`` (the Routh–Hurwitz test of ``Q(-z)``);
2. ``E(y) = |Q(iy)|^2 - |P(iy)|^2 ≥ 0`` for every real ``y`` (a sign count of the polynomial
   ``E`` on ``y > 0``, by its square-free factorisation and a Sturm sequence).

Both are exact for a rational `T`. For a floating-point `T`, a coefficient of ``Q`` counts as
zero by the zero count of [`stability_at_infinity`](@ref), a coefficient of ``E`` when its
magnitude is at most ``\\mathrm{rtol} \\cdot m^2``, with the same ``m``, and the decision is then
exact for the remaining coefficients. The default `rtol` is `16 * eps(T)` for a floating-point
`T` and `0` for a rational one. Returns a `Bool`.
"""
function isastable(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    P, Q, m = _zero_counted(stability_function(tab), rtol)
    E = _psub(_abs2_on_imaginary_axis(Q), _abs2_on_imaginary_axis(P))
    return _has_roots_in_right_half_plane_only(Q) && _isnonnegative(E, rtol * m^2)
end

"""
    islstable(tab::Tableau{T}; rtol) where {T <: Real}

Return whether `tab` is L-stable: [`isastable`](@ref) and ``R(∞) = 0``, that is
``\\deg P < \\deg Q`` by the zero count of [`stability_at_infinity`](@ref). The keyword `rtol`
is that of both. Returns a `Bool`.
"""
function islstable(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    P, Q, _ = _zero_counted(stability_function(tab), rtol)
    return isastable(tab; rtol) && _degree(P) < _degree(Q)
end

"""
    real_stability_interval(tab::Tableau{T}; rtol) where {T <: Real}

Return the supremum of the ``r ≥ 0`` with ``|R(x)| ≤ 1`` for every ``x ∈ [-r, 0]``, of type
`float(T)`, and `Inf` where no bound exists. A point where ``|R| = 1`` and ``|R| ≤ 1`` on both
sides does not end the interval.

It is decided from the coefficients of ``R = P/Q``: the interval ends at the first ``u > 0`` where
``Q(-u)^2 - P(-u)^2`` changes its sign to negative, found by a square-free factorisation and a
Sturm sequence in exact rational arithmetic, and an irrational endpoint is located by bisection to
the precision of `float(T)`. For a floating-point `T`, the coefficients of ``P`` and ``Q`` are
counted as zero as in [`stability_at_infinity`](@ref), and a coefficient of
``Q(-u)^2 - P(-u)^2`` when its magnitude is at most ``\\mathrm{rtol} \\cdot m^2``. The default
`rtol` is `16 * eps(T)` for a floating-point `T` and `0` for a rational one.
"""
function real_stability_interval(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    P, Q, m = _zero_counted(stability_function(tab), rtol)
    P, Q = _reflect(P), _reflect(Q)
    return _nonnegative_until(_psub(_pmul(Q, Q), _pmul(P, P)), rtol * m^2, float(T))
end

"""
    imaginary_stability_interval(tab::Tableau{T}; rtol) where {T <: Real}

Return the supremum of the ``r ≥ 0`` with ``|R(iy)| ≤ 1`` for every ``y ∈ [-r, r]``, of type
`float(T)`, and `Inf` where no bound exists. A point where ``|R| = 1`` and ``|R| ≤ 1`` on both
sides does not end the interval, and an interval where ``|R(iy)| > 1`` for every small
``y ≠ 0`` is exactly `0`.

It is decided from the coefficients of ``R = P/Q``, by the sign of
``E(y) = |Q(iy)|^2 - |P(iy)|^2`` on ``y > 0`` as in [`real_stability_interval`](@ref), with the
same zero count and the same default `rtol`.
"""
function imaginary_stability_interval(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <:
                                                                                       Real}
    P, Q, m = _zero_counted(stability_function(tab), rtol)
    E = _psub(_abs2_on_imaginary_axis(Q), _abs2_on_imaginary_axis(P))
    return _nonnegative_until(E, rtol * m^2, float(T))
end

# The coefficients, in ν, of the denominator D and of the numerators of det M and det(M + I),
# where M(ν) is the amplification matrix of one step on q̇ = p, ṗ = -q with the step ν, the
# tableau q acting on q̇ and the tableau p on ṗ. The stages solve (I + ν K) X = V, and
# M = I + ν W X, so that det(I + ν W X) = det(I + ν K + ν V W) / det(I + ν K).
function _amplification_coefficients(q::Tableau{T}, p::Tableau{T}) where {T}
    s = q.s
    o, z = ones(T, s), zeros(T, s)
    Z = zeros(T, s, s)
    K = [Z -Matrix{T}(q.a); Matrix{T}(p.a) Z]
    V = [o z; z o]
    W = [transpose(z) transpose(Vector{T}(q.b)); -transpose(Vector{T}(p.b)) transpose(z)]
    VW = V * W
    D = _det_coefficients(-K)
    N1 = _det_coefficients(-(K + VW))
    N2 = 4 .* _det_coefficients(-(K + VW / 2))
    return D, N1, N2
end

"""
    imaginary_stability_interval(tab::PartitionedTableau{T}; rtol) where {T <: Real}

Return the supremum of the ``r ≥ 0`` such that the amplification matrix ``M(ν)`` of one step on
the harmonic oscillator ``\\dot{q} = p``, ``\\dot{p} = -ω^2 q``, with ``ν = h ω``, the tableau
`tab.q` acting on ``\\dot{q}`` and the tableau `tab.p` on ``\\dot{p}``, has spectral radius at
most 1 for every ``ν ∈ [0, r]``; of type `float(T)`, and `Inf` where no bound exists. For a
tableau paired with itself, it equals the interval of that tableau.

The spectral radius of the real ``2 × 2`` matrix ``M`` is at most 1 exactly when
``|\\det M| ≤ 1``, ``\\det(I + M) ≥ 0`` and ``\\det(I - M) ≥ 0``. With ``\\det M = N_1/D`` and
``\\det(I + M) = N_2/D`` as quotients of polynomials in ``ν``, these are the polynomial
conditions ``D^2 - N_1^2 ≥ 0``, ``D N_2 ≥ 0`` and ``D (2D + 2N_1 - N_2) ≥ 0``, each decided as in
[`real_stability_interval`](@ref). For a floating-point `T`, a coefficient of ``D``, ``N_1`` or
``N_2`` counts as zero when its magnitude is at most ``\\mathrm{rtol} \\cdot m``, with ``m`` the
largest magnitude of their coefficients, and a coefficient of a condition when its magnitude is
at most ``\\mathrm{rtol} \\cdot m^2``. The default `rtol` is `16 * eps(T)` for a floating-point
`T` and `0` for a rational one.
"""
function imaginary_stability_interval(
        tab::PartitionedTableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    _check_finite(tab.q)
    _check_finite(tab.p)
    D, N1, N2 = _amplification_coefficients(tab.q, tab.p)
    m = max(maximum(abs, D), maximum(abs, N1), maximum(abs, N2))
    D, N1, N2 = _count_zeros(D, rtol * m), _count_zeros(N1, rtol * m),
    _count_zeros(N2, rtol * m)
    conditions = (_psub(_pmul(D, D), _pmul(N1, N1)), _pmul(D, N2),
        _pmul(D, _psub(_padd(2 .* D, 2 .* N1), N2)))
    return minimum(_nonnegative_until(g, rtol * m^2, float(T); touch = true)
    for g in conditions)
end

"""
    algebraic_stability_matrix(tab::Tableau{T})

Return the matrix ``M = B A + A^T B - b b^T``, ``B = \\mathrm{diag}(b)``, of `tab` as a
`Matrix{T}`, computed in the element type `T` with no tolerance. `tab` is algebraically stable
when ``b_i ≥ 0`` for every ``i`` and ``M`` is positive semidefinite
([`isalgebraicallystable`](@ref)); ``M = 0`` for a symplectic tableau.
"""
function algebraic_stability_matrix(tab::Tableau{T}) where {T}
    a, b = tab.a, tab.b
    return T[b[i] * a[i, j] + b[j] * a[j, i] - b[i] * b[j]
             for i in 1:(tab.s), j in 1:(tab.s)]
end

"""
    isalgebraicallystable(tab::Tableau{T}; rtol) where {T <: Real}

Return whether `tab` is algebraically stable: every weight ``b_i ≥ 0``, and the
[`algebraic_stability_matrix`](@ref) ``M`` is positive semidefinite. The semidefinite test is a
symmetric elimination that pivots on the largest diagonal entry; ``M`` is positive semidefinite
when each pivot is positive until the remaining block is zero.

The decision is exact for a rational `T`. For a floating-point `T`, a weight, a pivot or an entry
of the remaining block counts as zero when its magnitude is at most ``\\mathrm{rtol} \\cdot m``,
with ``m`` the largest magnitude of an entry of ``M`` and ``b``. The default `rtol` is
`16 * eps(T)` for a floating-point `T` and `0` for a rational one. Returns a `Bool`.
"""
function isalgebraicallystable(tab::Tableau{T}; rtol = _default_rtol(T)) where {T <: Real}
    _check_finite(tab)
    M = algebraic_stability_matrix(tab)
    tol = rtol * max(maximum(abs, M), maximum(abs, tab.b))
    all(>=(-tol), tab.b) || return false
    while !isempty(M)
        d, k = findmax([M[i, i] for i in axes(M, 1)])
        d <= tol && return all(x -> abs(x) <= tol, M)
        rest = [i for i in axes(M, 1) if i != k]
        M = M[rest, rest] - M[rest, k] * transpose(M[k, rest]) / d
    end
    return true
end

"""
    stability_region(tab::Tableau{T}; xlims = (-6.0, 2.0), ylims = (-4.0, 4.0), n = 201)

Return the data of the stability region of `tab` on a grid of the complex plane, as a
`NamedTuple` `(x, y, absR)`: `x` and `y` are the `n` equally spaced points of the closed
intervals `xlims` and `ylims`, and `absR[i, j] = |R(x[i] + im * y[j])|`, with `absR[i, j] = Inf`
wherever ``Q(x_i + i y_j) = 0``. All three hold `float(T)`. The stability region is the set where
`absR ≤ 1`.

Throws an `ArgumentError` when a limit is not finite, when a window is not increasing, or when
`n < 2`.
"""
function stability_region(tab::Tableau{T}; xlims = (-6.0, 2.0), ylims = (-4.0, 4.0),
        n = 201) where {T <: Real}
    all(isfinite, (xlims..., ylims...)) && xlims[1] < xlims[2] && ylims[1] < ylims[2] &&
    n >= 2 ||
        throw(ArgumentError("the window must have finite, increasing limits and n ≥ 2"))
    F = float(T)
    R = stability_function(tab)
    P, Q = F.(R.P), F.(R.Q)
    x = collect(range(F(xlims[1]), F(xlims[2]); length = n))
    y = collect(range(F(ylims[1]), F(ylims[2]); length = n))
    absR = [_abs_quotient(P, Q, complex(xi, yj)) for xi in x, yj in y]
    return (; x, y, absR)
end

function _abs_quotient(P, Q, z)
    q = evalpoly(z, Q)
    return iszero(q) ? oftype(abs(q), Inf) : abs(evalpoly(z, P) / q)
end
