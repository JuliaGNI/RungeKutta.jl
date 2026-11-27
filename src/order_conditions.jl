
# The default tolerance of the order-condition checks: exactly zero for an exact element type,
# so that the check compares exactly, and 16 machine epsilons of the element type otherwise.
_default_tolerance(::Type{T}) where {T} = 16eps(T)
_default_tolerance(::Type{T}) where {T <: Union{Integer, Rational}} = zero(T)

# The right-hand side x / k of an order condition. For an integer or rational x it is the exact
# rational x // k, because x / k converts an integer x to Float64, where integers above 2^53
# that differ are equal. It is formed in BigInt, so that isapprox subtracts it from the left-hand
# side without the OverflowError of checked Int64 rational arithmetic.
_divide(x, k) = x / k
_divide(x::Union{Integer, Rational}, k::Integer) = big(x) // k

function check_order_conditions_b(tab::Tableau{T}, k; atol = _default_tolerance(T),
        rtol = _default_tolerance(T)) where {T}
    b, c, s = tab.b, tab.c, tab.s
    isapprox(mapreduce(i -> b[i] * c[i]^(k-1), +, 1:s), _divide(one(T), k);
        atol = atol, rtol = rtol)
end

function satisfies_simplifying_assumption_b(tab::Tableau, σ = tab.s; kwargs...)
    all([check_order_conditions_b(tab, k; kwargs...) for k in 1:σ])
end

function check_order_conditions_c(tab::Tableau{T}, k; atol = _default_tolerance(T),
        rtol = _default_tolerance(T)) where {T}
    a, c, s = tab.a, tab.c, tab.s
    [isapprox(mapreduce(j -> a[i, j] * c[j]^(k-1), +, 1:s), _divide(c[i]^k, k); atol = atol, rtol = rtol)
     for i in 1:s]
end

function satisfies_simplifying_assumption_c(tab::Tableau, σ = tab.s; kwargs...)
    all(hcat([check_order_conditions_c(tab, k; kwargs...) for k in 1:σ]...))
end

function check_order_conditions_d(tab::Tableau{T}, k; atol = _default_tolerance(T),
        rtol = _default_tolerance(T)) where {T}
    a, b, c, s = tab.a, tab.b, tab.c, tab.s
    [isapprox(mapreduce(j -> b[j] * c[j]^(k-1) * a[j, i], +, 1:s),
         _divide(b[i] * (1 - c[i]^k), k); atol = atol, rtol = rtol) for i in 1:s]
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
