
@doc raw"""
The Radau IA coefficients are implicitly given by the so-called simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function radau_1_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau IA coefficients for one stage are not defined."))
    end
    solve_simplifying_assumption_d(
        radau_legendre_weights(T, s, Val(:left)), radau_legendre_nodes(T, s, Val(:left)))
end

radau_1_coefficients(s) = radau_1_coefficients(BigFloat, s)

@doc raw"""
The Radau IIA coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function radau_2_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Radau IIA coefficients for one stage are not defined."))
    end
    solve_simplifying_assumption_c(radau_legendre_nodes(T, s, Val(:right)))
end

radau_2_coefficients(s) = radau_2_coefficients(BigFloat, s)

function reference(::Val{:RadauIA})
    """
References:

    Byron Leonard Ehle
    On Padé approximations to the exponential function and a-stable methods for the numerical solution of initial value problems.
    Research Report CSRR 2010, Dept. AACS, University of Waterloo, 1969.
"""
end

@doc raw"""
Radau IA tableau with s stages

```julia
TableauRadauIA(::Type{T}, s)
TableauRadauIA(s) = TableauRadauIA(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

The nodes and weights are those of the **left** Radau-Legendre quadrature rule, i.e.
`QuadratureRules.radau_legendre_nodes(T, s, Val(:left))` and the corresponding weights: Radau
IA is the family that prescribes the *left* endpoint, so $c_1 = 0$. Equivalently the nodes are
the roots of $\frac{d^{s-1}}{dx^{s-1}} \big( x^s (x-1)^{s-1} \big)$. The coefficients follow
from the simplifying assumption $D(s)$, cf. [`radau_1_coefficients`](@ref).

Prescribing one endpoint costs one degree of exactness relative to Gauss, giving order $2s-1$.
Contrast [`TableauRadauIIA`](@ref), which prescribes the right endpoint instead.

""" * reference(Val(:RadauIA)) function TableauRadauIA(::Type{T}, s) where {T}
    Tableau{T}(:RadauIA, 2s-1, radau_1_coefficients(s),
        radau_legendre_weights(BigFloat, s, Val(:left)),
        radau_legendre_nodes(BigFloat, s, Val(:left)); R∞ = 0)
end

TableauRadauIA(s) = TableauRadauIA(Float64, s)

reference(::Val{:RadauIB}) = """
Reference:

    Sun Geng
    Construction of high order symplectic Runge-Kutta methods
    Journal of Computational Mathematics, Volume 11, Pages 250-260, 1993.
"""

"""
Radau IB tableau with s stages

```julia
TableauRadauIB(::Type{T}, s)
TableauRadauIB(s) = TableauRadauIB(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Coefficients are taken as ``a^B = \\frac{1}{2} ( a^A + \\bar{a}^A )`` where ``a^A`` are the coefficients
of the Radau IA method and ``\\bar{a}^A`` are computed such that the symplecticity conditions
``b_{i} \\bar{a}_{i,j} + \\bar{b}_{j} a_{j,i} = b_{i} \\bar{b}_{j}`` and ``b_{i} = \\bar{b}_i`` hold for
all ``1 \\le i,j \\le s``.

$(reference(Val(:RadauIB)))
"""
function TableauRadauIB(::Type{T}, s) where {T}
    a = radau_1_coefficients(BigFloat, s)
    b = radau_legendre_weights(BigFloat, s, Val(:left))
    ā = symplectic_conjugate_coefficients(a, b)

    Tableau{T}(:RadauIB, 2s-1, (a .+ ā) ./ 2, b,
        radau_legendre_nodes(BigFloat, s, Val(:left)); R∞ = 0)
end

TableauRadauIB(s) = TableauRadauIB(Float64, s)

function reference(::Val{:RadauIIA})
    """
References:

    Byron Leonard Ehle
    On Padé approximations to the exponential function and a-stable methods for the numerical solution of initial value problems.
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
end

@doc raw"""
Radau IIA tableau with s stages

```julia
TableauRadauIIA(::Type{T}, s)
TableauRadauIIA(s) = TableauRadauIIA(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

The nodes and weights are those of the **right** Radau-Legendre quadrature rule, i.e.
`QuadratureRules.radau_legendre_nodes(T, s, Val(:right))` and the corresponding weights: Radau
IIA is the family that prescribes the *right* endpoint, so $c_s = 1$. Equivalently the nodes
are the roots of $\frac{d^{s-1}}{dx^{s-1}} \big( x^{s-1} (x-1)^s \big)$. The coefficients
follow from the simplifying assumption $C(s)$, cf. [`radau_2_coefficients`](@ref).

Having the right endpoint among the nodes is what makes the method *stiffly accurate*, which is
why Radau IIA rather than [`TableauRadauIA`](@ref) is the workhorse for stiff and
differential-algebraic problems. Both have order $2s-1$.

""" * reference(Val(:RadauIIA)) function TableauRadauIIA(::Type{T}, s) where {T}
    Tableau{T}(:RadauIIA, 2s-1, radau_2_coefficients(s),
        radau_legendre_weights(BigFloat, s, Val(:right)),
        radau_legendre_nodes(BigFloat, s, Val(:right)); R∞ = 0)
end

TableauRadauIIA(s) = TableauRadauIIA(Float64, s)

reference(::Val{:RadauIIB}) = """
Reference:

    Sun Geng
    Construction of high order symplectic Runge-Kutta methods
    Journal of Computational Mathematics, Volume 11, Pages 250-260, 1993.
"""

"""
Radau IIB tableau with s stages

```julia
TableauRadauIIB(::Type{T}, s)
TableauRadauIIB(s) = TableauRadauIIB(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Coefficients are taken as ``a^B = \\frac{1}{2} ( a^A + \\bar{a}^A )`` where ``a^A`` are the coefficients
of the Radau IIA method and ``\\bar{a}^AV are computed such that the symplecticity conditions
``b_{i} \\bar{a}_{i,j} + \\bar{b}_{j} a_{j,i} = b_{i} \\bar{b}_{j}`` and ``b_{i} = \\bar{b}_i`` hold for
all ``1 \\le i,j \\le s``.

$(reference(Val(:RadauIIB)))
"""
function TableauRadauIIB(::Type{T}, s) where {T}
    a = radau_2_coefficients(BigFloat, s)
    b = radau_legendre_weights(BigFloat, s, Val(:right))
    ā = symplectic_conjugate_coefficients(a, b)

    Tableau{T}(:RadauIIB, 2s-1, (a .+ ā) ./ 2, b,
        radau_legendre_nodes(BigFloat, s, Val(:right)); R∞ = 0)
end

TableauRadauIIB(s) = TableauRadauIIB(Float64, s)
