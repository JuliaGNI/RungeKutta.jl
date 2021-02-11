
@doc raw"""
The s-stage Radau IA nodes are defined as the roots of the following polynomial of degree $s$:
```math
\frac{d^{s-1}}{dx^{s-1}} \big( x^s (x - 1)^{s-1} \big) .
```
"""
function get_radau_1_nodes(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau nodes for one stage are not defined."))
    end

    D(k) = Polynomials.derivative(Polynomial(T[0, 1])^k * Polynomial(T[-1, 1])^(k-1), k-1)
    c = sort(T.(Polynomials.roots(D(s))))
    c[begin] = 0; c
end

get_radau_1_nodes(s) = get_radau_1_nodes(BigFloat, s)


@doc raw"""
The s-stage Radau IIA nodes are defined as the roots of the following polynomial of degree $s$:
```math
\frac{d^{s-1}}{dx^{s-1}} \big( x^{s-1} (x - 1)^s \big) .
```
"""
function get_radau_2_nodes(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau nodes for one stage are not defined."))
    end

    D(k) = Polynomials.derivative(Polynomial(T[0, 1])^(k-1) * Polynomial(T[-1, 1])^k, k-1)
    c = sort(T.(Polynomials.roots(D(s))))
    c[end] = 1; c
end

get_radau_2_nodes(s) = get_radau_2_nodes(BigFloat, s)


@doc raw"""
The Radau IA weights are implicitly given by the so-called simplifying assumption $B(s)$:
```math
\sum \limits_{j=1}^{s} b_{j} c_{j}^{k-1} = \frac{1}{k}  \qquad k = 1 , \, ... , \, s .
```
"""
function get_radau_1_weights(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau weights for one stage are not defined."))
    end
    solve_simplifying_assumption_b(get_radau_1_nodes(T,s))
end

get_radau_1_weights(s) = get_radau_1_weights(BigFloat, s)


@doc raw"""
The Radau IIA weights are implicitly given by the so-called simplifying assumption $B(s)$:
```math
\sum \limits_{j=1}^{s} b_{j} c_{j}^{k-1} = \frac{1}{k}  \qquad k = 1 , \, ... , \, s .
```
"""
function get_radau_2_weights(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau weights for one stage are not defined."))
    end
    solve_simplifying_assumption_b(get_radau_2_nodes(T,s))
end

get_radau_2_weights(s) = get_radau_2_weights(BigFloat, s)


@doc raw"""
The Radau IA coefficients are implicitly given by the so-called simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_radau_1_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau IIA coefficients for one stage are not defined."))
    end
    solve_simplifying_assumption_d(get_radau_1_weights(T,s), get_radau_1_nodes(T,s))
end

get_radau_1_coefficients(s) = get_radau_1_coefficients(BigFloat, s)


@doc raw"""
The Radau IIA coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_radau_2_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau IIA coefficients for one stage are not defined."))
    end
    solve_simplifying_assumption_c(get_radau_2_nodes(T,s))
end

get_radau_2_coefficients(s) = get_radau_2_coefficients(BigFloat, s)


"""
Radau IA tableau with s stages

```julia
TableauRadauIA(::Type{T}, s)
TableauRadauIA(s) = TableauRadauIA(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

References:

    Byron Leonard Ehle
    On Padé approximations to the exponential function and a-stable methods for the numerical solution of initial value problems.
    Research Report CSRR 2010, Dept. AACS, University of Waterloo, 1969.

"""
function TableauRadauIA(::Type{T}, s) where {T}
    Tableau{T}(Symbol("RadauIA($s)"), 2s-1, get_radau_1_coefficients(s), get_radau_1_weights(s), get_radau_1_nodes(s))
end

TableauRadauIA(s) = TableauRadauIA(Float64, s)


@doc raw"""
Radau IB tableau with s stages

```julia
TableauRadauIB(::Type{T}, s)
TableauRadauIB(s) = TableauRadauIB(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Coefficients are taken as $a^B = \frac{1}{2} ( a^A + \bar{a}^A )$ where $a^A$ are the coefficients
of the Radau IA method and $\bar{a}^A$ are computed such that the symplecticity conditions
$b_{i} \bar{a}_{i,j} + \bar{b}_{j} a_{j,i} = b_{i} \bar{b}_{j}$ and $b_{i} = \bar{b}_i$ hold for
all $1 \le i,j \le s$.

Reference:

    Sun Geng
    Construction of high order symplectic Runge-Kutta methods
    Journal of Computational Mathematics, Volume 11, Pages 250-260, 1993.

"""
function TableauRadauIB(::Type{T}, s) where {T}
    a = get_radau_1_coefficients(BigFloat,s)
    b = get_radau_1_weights(BigFloat,s)
    ā = get_symplectic_conjugate_coefficients(a,b)

    Tableau{T}(Symbol("RadauIB($s)"), 2s-1, (a .+ ā) ./ 2, b, get_radau_1_nodes(s))
end

TableauRadauIB(s) = TableauRadauIB(Float64, s)


"""
Radau IIA tableau with s stages

```julia
TableauRadauIIA(::Type{T}, s)
TableauRadauIIA(s) = TableauRadauIIA(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

References:

    Byron Leonard Ehle
    On Padé approximations to the exponential function and a-stable methods for the numerical solution of initial value problems.
    Research Report CSRR 2010, Dept. AACS, University of Waterloo, 1969.

    Owe Axelsson.
    A class of A-stable methods.
    BIT, Volume 9, Pages 185-199, 1969.
    doi: 10.1007/BF01946812.

    Ernst Hairer and Gerhard Wanner.
    Radau Methods.
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_139.

"""
function TableauRadauIIA(::Type{T}, s) where {T}
    Tableau{T}(Symbol("RadauIIA($s)"), 2s-1, get_radau_2_coefficients(s), get_radau_2_weights(s), get_radau_2_nodes(s))
end

TableauRadauIIA(s) = TableauRadauIIA(Float64, s)


@doc raw"""
Radau IIB tableau with s stages

```julia
TableauRadauIIB(::Type{T}, s)
TableauRadauIIB(s) = TableauRadauIIB(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Coefficients are taken as $a^B = \frac{1}{2} ( a^A + \bar{a}^A )$ where $a^A$ are the coefficients
of the Radau IIA method and $\bar{a}^A$ are computed such that the symplecticity conditions
$b_{i} \bar{a}_{i,j} + \bar{b}_{j} a_{j,i} = b_{i} \bar{b}_{j}$ and $b_{i} = \bar{b}_i$ hold for
all $1 \le i,j \le s$.

Reference:

    Sun Geng
    Construction of high order symplectic Runge-Kutta methods
    Journal of Computational Mathematics, Volume 11, Pages 250-260, 1993.

"""
function TableauRadauIIB(::Type{T}, s) where {T}
    a = get_radau_2_coefficients(BigFloat,s)
    b = get_radau_2_weights(BigFloat,s)
    ā = get_symplectic_conjugate_coefficients(a,b)

    Tableau{T}(Symbol("RadauIIB($s)"), 2s-1, (a .+ ā) ./ 2, b, get_radau_2_nodes(s))
end

TableauRadauIIB(s) = TableauRadauIIB(Float64, s)
