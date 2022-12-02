
import CompactBasisFunctions: Lagrange
import LinearAlgebra


@doc raw"""
The s-stage Lobatto nodes are defined as the roots of the following polynomial of degree $s$:
```math
\frac{d^{s-2}}{dx^{s-2}} \big( (x - x^2)^{s-1} \big) .
```
"""
function get_lobatto_nodes(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto nodes for one stage are not defined."))
    end

    D(k) = Polynomials.derivative(Polynomial(T[0, 1, -1])^(k-1), k-2)
    c = sort(T.(Polynomials.roots(D(s))))
    c[begin] = 0; c[end] = 1; c
end

get_lobatto_nodes(s) = get_lobatto_nodes(BigFloat, s)


@doc raw"""
The Lobatto weights can be explicitly computed by the formula
```math
b_j = \frac{1}{s (s-1) P_{s-1}(2 c_j - 1)^2} \qquad j = 1 , \, ... , \, s ,
```
where $P_k$ is the $k$th Legendre polynomial, given by
```math
P_k (x) = \frac{1}{k! 2^k} \big( \frac{d^k}{dx^k} (x^2 - 1)^k \big) .
```
"""
function get_lobatto_weights(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto weights for one stage are not defined."))
    end

    P(k,x) = Polynomials.derivative(Polynomial(T[-1, 0, 1])^k, k)(x) / LinearAlgebra.factorial(k) / 2^k
    c = get_lobatto_nodes(T,s)
    b = [ 1 / ( s*(s-1) * P(s-1, 2c[i] - 1)^2 ) for i in 1:s ]
end

get_lobatto_weights(s) = get_lobatto_weights(BigFloat, s)

"""
```julia
get_lobatto_nullvector(::Type, s; normalize=false)
get_lobatto_nullvector(s; kwargs...)
```

Computes the nullvector of the matrix containing the derivatives of the
Lagrange basis on the `s` Lobatto nodes evaluated on these nodes.
"""
function get_lobatto_nullvector(::Type{T}, s; normalize=false) where {T}
    if s == 1
        throw(ErrorException("Lobatto nullvector for one stage is not defined."))
    end

    q = get_lobatto_nodes(s)
    l = Lagrange(q)
    v = [l'[x, j] for x in q, j in eachindex(l)]
    w = LinearAlgebra.nullspace(v')[:,begin]
    normalize ? T.(LinearAlgebra.normalize(w) .* sign(w[begin])) : T.(w)
end

get_lobatto_nullvector(s; kwargs...) = get_lobatto_nullvector(Float64, s; kwargs...)


@doc raw"""
The Lobatto IIIA coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_lobatto_a_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIA coefficients for one stage are not defined."))
    end
    solve_simplifying_assumption_c(get_lobatto_nodes(T,s))
end

@doc raw"""
The Lobatto IIIB coefficients are implicitly given by the so-called simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_lobatto_b_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIB coefficients for one stage are not defined."))
    end

    solve_simplifying_assumption_d(get_lobatto_weights(T,s), get_lobatto_nodes(T,s))
end

@doc raw"""
The Lobatto IIIC coefficients are determined by setting $a_{i,1} = b_1$ and
solving the so-called simplifying assumption $C(s-1)$, given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s-1 ,
```
for $a_{i,j}$ with $i = 1, ..., s$ and $j = 2, ..., s$.
"""
function get_lobatto_c_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIC coefficients for one stage are not defined."))
    end

    b = get_lobatto_weights(T,s)
    c = get_lobatto_nodes(T,s)
    M = [ c[j]^(k-1) for k in 1:s-1, j in 2:s ]
    
    row(i) = begin
        r = [ c[i]^k / T(k) - c[1]^(k-1) * b[1] for k in 1:s-1 ]
        M \ r
    end
    
    hcat(b[1] * ones(T,s), vcat([transpose(row(i)) for i in 1:s]...))
end

@doc raw"""
The Lobatto IIIC̄ coefficients are determined by setting $a_{i,s} = 0$ and
solving the so-called simplifying assumption $C(s-1)$, given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s-1 ,
```
for $a_{i,j}$ with $i = 1, ..., s$ and $j = 1, ..., s-1$.
"""
function get_lobatto_c̄_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIC̄ coefficients for one stage are not defined."))
    end

    c = get_lobatto_nodes(T,s)
    M = [ c[j]^(k-1) for k in 1:s-1, j in 1:s-1 ]
    
    row(i) = begin
        r = [ c[i]^k / T(k) for k in 1:s-1 ]
        M \ r
    end
    
    hcat(vcat([transpose(row(i)) for i in 1:s]...), zeros(T,s))
end

get_lobatto_d_coefficients(::Type{T}, s) where {T} = (get_lobatto_c_coefficients(T,s) .+ get_lobatto_c̄_coefficients(T,s)) ./ 2

get_lobatto_e_coefficients(::Type{T}, s) where {T} = (get_lobatto_a_coefficients(T,s) .+ get_lobatto_b_coefficients(T,s)) ./ 2

function get_lobatto_f_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIF coefficients for one stage are not defined."))
    end

    c = get_lobatto_nodes(T,s)
    M = [ 1 / T(k + j - 1) for k in 1:s, j in 1:s ]
    r = [ 1 / T(s) / T(s + k) for k in 1:s ]
    α = M \ r
    
    Vₛ = [ c[i]^(j-1) for i in 1:s, j in 1:s ]
    Aₛ = zeros(T, s, s)
    for i in 2:s
       Aₛ[i,i-1] = 1 / T(i-1)
    end
    Aₛ[:,s] = α
    
    Vₛ * Aₛ * inv(Vₛ)
end

function get_lobatto_g_coefficients(::Type{T}, s) where {T}
    a = get_lobatto_f_coefficients(T,s)
    b = get_lobatto_weights(T,s)
    ā = symplectic_conjugate_coefficients(a,b)
    return (a .+ ā) ./ 2
end


get_lobatto_a_coefficients(s) = get_lobatto_a_coefficients(BigFloat, s)
get_lobatto_b_coefficients(s) = get_lobatto_b_coefficients(BigFloat, s)
get_lobatto_c_coefficients(s) = get_lobatto_c_coefficients(BigFloat, s)
get_lobatto_c̄_coefficients(s) = get_lobatto_c̄_coefficients(BigFloat, s)
get_lobatto_d_coefficients(s) = get_lobatto_d_coefficients(BigFloat, s)
get_lobatto_e_coefficients(s) = get_lobatto_e_coefficients(BigFloat, s)
get_lobatto_f_coefficients(s) = get_lobatto_f_coefficients(BigFloat, s)
get_lobatto_g_coefficients(s) = get_lobatto_g_coefficients(BigFloat, s)


RungeKutta.reference(::Val{:LobattoIII}) = """
References:

    John C. Butcher.
    Integration processes based on Radau quadrature formulas
    Mathematics of Computation, Volume 18, Pages 233-244, 1964.
    doi: 10.1090/S0025-5718-1964-0165693-1.

    Laurent O. Jay.
    Lobatto Methods.
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_123.
"""

"""
Lobatto III tableau with s stages

```julia
TableauLobattoIII(::Type{T}, s)
TableauLobattoIII(s) = TableauLobattoIII(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Sometimes this tableau is also referred to as Lobatto IIIC*.

$(reference(Val(:LobattoIII)))
"""
function TableauLobattoIII(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIII, 2s-2, get_lobatto_c̄_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^(s+1))
end


RungeKutta.reference(::Val{:LobattoIIIA}) = """
References:

    Byron Leonard Ehle
    On Padé approximations to the exponential function and a-stable methods for the numerical solution of initial value problems.
    Research Report CSRR 2010, Dept. AACS, University of Waterloo, 1969.

    Laurent O. Jay.
    Lobatto Methods.
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_123
"""

"""
Lobatto IIIA tableau with s stages

```julia
TableauLobattoIIIA(::Type{T}, s)
TableauLobattoIIIA(s) = TableauLobattoIIIA(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

$(reference(Val(:LobattoIIIA)))
"""
function TableauLobattoIIIA(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIIIA, 2s-2, get_lobatto_a_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^(s+1))
end


"""
Lobatto IIIĀ tableau with s stages

```julia
TableauLobattoIIIĀ(::Type{T}, s)
TableauLobattoIIIĀ(s) = TableauLobattoIIIĀ(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Lobatto IIIĀ tableau is the conjugate symplectic to [`TableauLobattoIIIA`](@ref).
On paper, its coefficients are identical to [`TableauLobattoIIIB`](@ref), however, they are computed
by the symplecticity condition and not by the formula for Lobatto IIIB and thus the numerical
values are slightly different.
"""
function TableauLobattoIIIĀ(::Type{T}, s) where {T}
    a = get_lobatto_a_coefficients(s)
    b = get_lobatto_weights(s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIĀ, 2s-2, ā, b, get_lobatto_nodes(s); R∞=(-1)^(s+1))
end


RungeKutta.reference(::Val{:LobattoIIIB}) = """
References:

    Byron Leonard Ehle.
    On Padé approximations to the exponential function and a-stable methods for the numerical solution of initial value problems.
    Research Report CSRR 2010, Dept. AACS, University of Waterloo, 1969.

    Laurent O. Jay.
    Lobatto Methods.
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_123.
"""

"""
Lobatto IIIB tableau with s stages

```julia
TableauLobattoIIIB(::Type{T}, s)
TableauLobattoIIIB(s) = TableauLobattoIIIB(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

$(reference(Val(:LobattoIIIB)))
"""
function TableauLobattoIIIB(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIIIB, 2s-2, get_lobatto_b_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^(s+1))
end


"""
Lobatto IIIB̄ tableau with s stages

```julia
TableauLobattoIIIB̄(::Type{T}, s)
TableauLobattoIIIB̄(s) = TableauLobattoIIIB̄(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Lobatto IIIB̄ tableau is the conjugate symplectic to [`TableauLobattoIIIB`](@ref).
On paper, its coefficients are identical to [`TableauLobattoIIIA`](@ref), however, they are computed
by the symplecticity condition and not by the formula for Lobatto IIIA and thus the numerical
values are slightly different.
"""
function TableauLobattoIIIB̄(::Type{T}, s) where {T}
    a = get_lobatto_b_coefficients(s)
    b = get_lobatto_weights(s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIB̄, 2s-2, ā, b, get_lobatto_nodes(s); R∞=(-1)^(s+1))
end


RungeKutta.reference(::Val{:LobattoIIIC}) = """
References:

    F. H. Chipman.
    A-stable Runge-Kutta processes.
    BIT, Volume 11, Pages 384-388, 1971.
    doi: 10.1007/BF01939406.

    Laurent O. Jay.
    Lobatto Methods.
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_123.
"""

"""
Lobatto IIIC tableau with s stages

```julia
TableauLobattoIIIC(::Type{T}, s)
TableauLobattoIIIC(s) = TableauLobattoIIIC(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

$(reference(Val(:LobattoIIIC)))
"""
function TableauLobattoIIIC(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIIIC, 2s-2, get_lobatto_c_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^(s+1))
end


"""
Lobatto IIIC̄ tableau with s stages

```julia
TableauLobattoIIIC̄(::Type{T}, s)
TableauLobattoIIIC̄(s) = TableauLobattoIIIC̄(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Lobatto IIIC̄ tableau is the conjugate symplectic to [`TableauLobattoIIIC`](@ref).
On paper, its coefficients are identical to [`TableauLobattoIII`](@ref), however, they are computed
by the symplecticity condition and not by the formula for Lobatto III and thus the numerical
values are slightly different.
"""
function TableauLobattoIIIC̄(::Type{T}, s) where {T}
    a = get_lobatto_c_coefficients(s)
    b = get_lobatto_weights(s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIC̄, 2s-2, ā, b, get_lobatto_nodes(s); R∞=(-1)^(s+1))
end


RungeKutta.reference(::Val{:LobattoIIID}) = """
References:

    R.P.K. Chan.
    On symmetric Runge-Kutta methods of high order.
    Computing, Volume 45, Pages 301-309, 1990.
    doi: 10.1007/BF02238798

    Laurent O. Jay.
    Lobatto Methods.
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_123.
"""

"""
Lobatto IIID tableau with s stages

```julia
TableauLobattoIIID(::Type{T}, s)
TableauLobattoIIID(s) = TableauLobattoIIID(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

$(reference(Val(:LobattoIIID)))
"""
function TableauLobattoIIID(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIIID, 2s-2, get_lobatto_d_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^s)
end


"""
Lobatto IIID̄ tableau with s stages

```julia
TableauLobattoIIID̄(::Type{T}, s)
TableauLobattoIIID̄(s) = TableauLobattoIIID̄(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Lobatto IIID̄ tableau is the conjugate symplectic to [`TableauLobattoIIID`](@ref).
On paper, the coefficients of the Lobatto IIID tableau are symplectic, however, the Lobatto IIID̄
coefficients are computed by the symplecticity condition and not by the formula for Lobatto IIID
and thus the numerical values are slightly different.
"""
function TableauLobattoIIID̄(::Type{T}, s) where {T}
    a = get_lobatto_d_coefficients(s)
    b = get_lobatto_weights(s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIID̄, 2s-2, ā, b, get_lobatto_nodes(s); R∞=(-1)^s)
end


RungeKutta.reference(::Val{:LobattoIIIE}) = """
References:

    R.P.K. Chan.
    On symmetric Runge-Kutta methods of high order.
    Computing, Volume 45, Pages 301-309, 1990.
    doi: 10.1007/BF02238798

    Laurent O. Jay.
    Lobatto Methods.
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_123.
"""

"""
Lobatto IIIE tableau with s stages

```julia
TableauLobattoIIIE(::Type{T}, s)
TableauLobattoIIIE(s) = TableauLobattoIIIE(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

$(reference(Val(:LobattoIIIE)))
"""
function TableauLobattoIIIE(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIIIE, 2s-2, get_lobatto_e_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^s)
end


"""
Lobatto IIIĒ tableau with s stages

```julia
TableauLobattoIIIĒ(::Type{T}, s)
TableauLobattoIIIĒ(s) = TableauLobattoIIIĒ(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Lobatto IIIĒ tableau is the conjugate symplectic to [`TableauLobattoIIIE`](@ref).
On paper, the coefficients of the Lobatto IIIE tableau are symplectic, however, the Lobatto IIIĒ
coefficients are computed by the symplecticity condition and not by the formula for Lobatto IIIE
and thus the numerical values are slightly different.
"""
function TableauLobattoIIIĒ(::Type{T}, s) where {T}
    a = get_lobatto_e_coefficients(s)
    b = get_lobatto_weights(s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIĒ, 2s-2, ā, b, get_lobatto_nodes(s); R∞=(-1)^s)
end


RungeKutta.reference(::Val{:LobattoIIIF}) = """
References:

    Wang Fangzong and Liao Xiaobing.
    A Class of Lobatto Methods of Order 2s.
    Journal of Applied Mathematics, Volume 46, Pages 6-10, 2016.
"""

"""
Lobatto IIIF tableau with s stages

```julia
TableauLobattoIIIF(::Type{T}, s)
TableauLobattoIIIF(s) = TableauLobattoIIIF(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

$(reference(Val(:LobattoIIIF)))
"""
function TableauLobattoIIIF(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIIIF, 2s,   get_lobatto_f_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^s)
end


"""
Lobatto IIIF̄ tableau with s stages

```julia
TableauLobattoIIIF̄(::Type{T}, s)
TableauLobattoIIIF̄(s) = TableauLobattoIIIF̄(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

The Lobatto IIIF̄ tableau is the conjugate symplectic to [`TableauLobattoIIIF`](@ref).
"""
function TableauLobattoIIIF̄(::Type{T}, s) where {T}
    a = get_lobatto_f_coefficients(s)
    b = get_lobatto_weights(s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIF̄, 2s, ā, b, get_lobatto_nodes(s); R∞=(-1)^s)
end


@doc raw"""
Lobatto IIIG tableau with s stages

```julia
TableauLobattoIIIG(::Type{T}, s)
TableauLobattoIIIG(s) = TableauLobattoIIIG(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Symplectizied algorithm for [`TableauLobattoIIIF`](@ref)

Coefficients are taken as $a^G = \frac{1}{2} ( a^F + \bar{a}^F )$ where the coefficients $\bar{a}^F$ are computed such that
the symplecticity conditions $b_{i} \bar{a}_{i,j} + \bar{b}_{j} a_{j,i} = b_{i} \bar{b}_{j}$ and $b_{i} = \bar{b}_i$ hold
for all $1 \le i,j \le s$.
"""
function TableauLobattoIIIG(::Type{T}, s) where {T}
    Tableau{T}(:LobattoIIIG, 2s,   get_lobatto_g_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s); R∞=(-1)^s)
end


"""
Lobatto IIIḠ tableau with s stages

```julia
TableauLobattoIIIḠ(::Type{T}, s)
TableauLobattoIIIḠ(s) = TableauLobattoIIIḠ(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

Lobatto IIIḠ tableau is the conjugate symplectic to [`TableauLobattoIIIG`](@ref).
On paper, the coefficients of the Lobatto IIIG tableau are symplectic, however, the Lobatto IIIḠ
coefficients are computed by the symplecticity condition and not by the formula for Lobatto IIIG
and thus the numerical values are slightly different.
"""
function TableauLobattoIIIḠ(::Type{T}, s) where {T}
    a = get_lobatto_g_coefficients(s)
    b = get_lobatto_weights(s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIḠ, 2s, ā, b, get_lobatto_nodes(s); R∞=(-1)^s)
end


TableauLobattoIII(s) = TableauLobattoIII(Float64, s)
TableauLobattoIIIA(s) = TableauLobattoIIIA(Float64, s)
TableauLobattoIIIĀ(s) = TableauLobattoIIIĀ(Float64, s)
TableauLobattoIIIB(s) = TableauLobattoIIIB(Float64, s)
TableauLobattoIIIB̄(s) = TableauLobattoIIIB̄(Float64, s)
TableauLobattoIIIC(s) = TableauLobattoIIIC(Float64, s)
TableauLobattoIIIC̄(s) = TableauLobattoIIIC̄(Float64, s)
TableauLobattoIIID(s) = TableauLobattoIIID(Float64, s)
TableauLobattoIIID̄(s) = TableauLobattoIIID̄(Float64, s)
TableauLobattoIIIE(s) = TableauLobattoIIIE(Float64, s)
TableauLobattoIIIĒ(s) = TableauLobattoIIIĒ(Float64, s)
TableauLobattoIIIF(s) = TableauLobattoIIIF(Float64, s)
TableauLobattoIIIF̄(s) = TableauLobattoIIIF̄(Float64, s)
TableauLobattoIIIG(s) = TableauLobattoIIIG(Float64, s)
TableauLobattoIIIḠ(s) = TableauLobattoIIIḠ(Float64, s)
