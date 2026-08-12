
import CompactBasisFunctions: Lagrange


"""
```julia
lobatto_nullvector(::Type, s)
lobatto_nullvector(s)
```

Computes the nullvector of the matrix containing the derivatives of the
Lagrange basis on the `s` Lobatto nodes evaluated on these nodes.

The result is of unit length with a positive first entry, cf. [`_nullvector`](@ref),
so it is determined by the Lobatto nodes alone and not by the factorisation used to
obtain it.
"""
function lobatto_nullvector(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto nullvector for one stage is not defined."))
    end

    q = lobatto_legendre_nodes(BigFloat, s)
    l = Lagrange(q)
    v = [l'[x, j] for x in q, j in eachindex(l)]
    T.(_nullvector(v'))
end

lobatto_nullvector(s) = lobatto_nullvector(Float64, s)


@doc raw"""
The Lobatto IIIA coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function lobatto_a_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIA coefficients for one stage are not defined."))
    end
    solve_simplifying_assumption_c(lobatto_legendre_nodes(T, s))
end

@doc raw"""
The Lobatto IIIB coefficients are implicitly given by the so-called simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function lobatto_b_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIB coefficients for one stage are not defined."))
    end

    solve_simplifying_assumption_d(lobatto_legendre_weights(T, s), lobatto_legendre_nodes(T, s))
end

@doc raw"""
The Lobatto IIIC coefficients are determined by setting $a_{i,1} = b_1$ and
solving the so-called simplifying assumption $C(s-1)$, given by
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s-1 ,
```
for $a_{i,j}$ with $i = 1, ..., s$ and $j = 2, ..., s$.
"""
function lobatto_c_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIC coefficients for one stage are not defined."))
    end

    b = lobatto_legendre_weights(T, s)
    c = lobatto_legendre_nodes(T, s)
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
function lobatto_c̄_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIC̄ coefficients for one stage are not defined."))
    end

    c = lobatto_legendre_nodes(T, s)
    M = [ c[j]^(k-1) for k in 1:s-1, j in 1:s-1 ]
    
    row(i) = begin
        r = [ c[i]^k / T(k) for k in 1:s-1 ]
        M \ r
    end
    
    hcat(vcat([transpose(row(i)) for i in 1:s]...), zeros(T,s))
end

lobatto_d_coefficients(::Type{T}, s) where {T} = (lobatto_c_coefficients(T,s) .+ lobatto_c̄_coefficients(T,s)) ./ 2

lobatto_e_coefficients(::Type{T}, s) where {T} = (lobatto_a_coefficients(T,s) .+ lobatto_b_coefficients(T,s)) ./ 2

function lobatto_f_coefficients(::Type{T}, s) where {T}
    if s == 1
        throw(ErrorException("Lobatto IIIF coefficients for one stage are not defined."))
    end

    c = lobatto_legendre_nodes(T, s)
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

function lobatto_g_coefficients(::Type{T}, s) where {T}
    a = lobatto_f_coefficients(T,s)
    b = lobatto_legendre_weights(T, s)
    ā = symplectic_conjugate_coefficients(a,b)
    return (a .+ ā) ./ 2
end


lobatto_a_coefficients(s) = lobatto_a_coefficients(BigFloat, s)
lobatto_b_coefficients(s) = lobatto_b_coefficients(BigFloat, s)
lobatto_c_coefficients(s) = lobatto_c_coefficients(BigFloat, s)
lobatto_c̄_coefficients(s) = lobatto_c̄_coefficients(BigFloat, s)
lobatto_d_coefficients(s) = lobatto_d_coefficients(BigFloat, s)
lobatto_e_coefficients(s) = lobatto_e_coefficients(BigFloat, s)
lobatto_f_coefficients(s) = lobatto_f_coefficients(BigFloat, s)
lobatto_g_coefficients(s) = lobatto_g_coefficients(BigFloat, s)


reference(::Val{:LobattoIII}) = """
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
    Tableau{T}(:LobattoIII, 2s-2, lobatto_c̄_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^(s+1))
end


reference(::Val{:LobattoIIIA}) = """
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
    Tableau{T}(:LobattoIIIA, 2s-2, lobatto_a_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^(s+1))
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
    a = lobatto_a_coefficients(s)
    b = lobatto_legendre_weights(BigFloat, s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIĀ, 2s-2, ā, b, lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^(s+1))
end


reference(::Val{:LobattoIIIB}) = """
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
    Tableau{T}(:LobattoIIIB, 2s-2, lobatto_b_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^(s+1))
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
    a = lobatto_b_coefficients(s)
    b = lobatto_legendre_weights(BigFloat, s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIB̄, 2s-2, ā, b, lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^(s+1))
end


reference(::Val{:LobattoIIIC}) = """
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
    Tableau{T}(:LobattoIIIC, 2s-2, lobatto_c_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^(s+1))
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
    a = lobatto_c_coefficients(s)
    b = lobatto_legendre_weights(BigFloat, s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIC̄, 2s-2, ā, b, lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^(s+1))
end


reference(::Val{:LobattoIIID}) = """
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
    Tableau{T}(:LobattoIIID, 2s-2, lobatto_d_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
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
    a = lobatto_d_coefficients(s)
    b = lobatto_legendre_weights(BigFloat, s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIID̄, 2s-2, ā, b, lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
end


reference(::Val{:LobattoIIIE}) = """
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
    Tableau{T}(:LobattoIIIE, 2s-2, lobatto_e_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
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
    a = lobatto_e_coefficients(s)
    b = lobatto_legendre_weights(BigFloat, s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIĒ, 2s-2, ā, b, lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
end


reference(::Val{:LobattoIIIF}) = """
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
    Tableau{T}(:LobattoIIIF, 2s,   lobatto_f_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
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
    a = lobatto_f_coefficients(s)
    b = lobatto_legendre_weights(BigFloat, s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIF̄, 2s, ā, b, lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
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
    Tableau{T}(:LobattoIIIG, 2s,   lobatto_g_coefficients(s), lobatto_legendre_weights(BigFloat, s), lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
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
    a = lobatto_g_coefficients(s)
    b = lobatto_legendre_weights(BigFloat, s)
    ā = symplectic_conjugate_coefficients(a,b)
    Tableau{T}(:LobattoIIIḠ, 2s, ā, b, lobatto_legendre_nodes(BigFloat, s); R∞=(-1)^s)
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
