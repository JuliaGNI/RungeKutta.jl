
import LinearAlgebra: factorial
import Polynomials
import Polynomials: Polynomial


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

    P(k,x) = Polynomials.derivative(Polynomial(T[-1, 0, 1])^k, k)(x) / factorial(k) / 2^k
    c = get_lobatto_nodes(T,s)
    b = [ 1 / ( s*(s-1) * P(s-1, 2c[i] - 1)^2 ) for i in 1:s ]
end

get_lobatto_weights(s) = get_lobatto_weights(BigFloat, s)


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
    ā = get_symplectic_conjugate_coefficients(a,b)
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


"Lobatto IIIA tableau with s stages"
function TableauLobattoIIIA(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIIA($s)"), 2s-2, get_lobatto_a_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIB tableau with s stages"
function TableauLobattoIIIB(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIIB($s)"), 2s-2, get_lobatto_b_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIC tableau with s stages"
function TableauLobattoIIIC(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIIC($s)"), 2s-2, get_lobatto_c_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIC̄ tableau with s stages"
function TableauLobattoIIIC̄(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIIC̄($s)"), 2s-2, get_lobatto_c̄_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIID tableau with s stages"
function TableauLobattoIIID(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIID($s)"), 2s-2, get_lobatto_d_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIE tableau with s stages"
function TableauLobattoIIIE(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIIE($s)"), 2s-2, get_lobatto_e_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIF tableau with s stages"
function TableauLobattoIIIF(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIIF($s)"), 2s,   get_lobatto_f_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end

"Lobatto IIIG tableau with s stages"
function TableauLobattoIIIG(::Type{T}, s) where {T}
    Tableau{T}(Symbol("LobattoIIIG($s)"), 2s,   get_lobatto_g_coefficients(s), get_lobatto_weights(s), get_lobatto_nodes(s))
end


TableauLobattoIIIA(s) = TableauLobattoIIIA(Float64, s)
TableauLobattoIIIB(s) = TableauLobattoIIIB(Float64, s)
TableauLobattoIIIC(s) = TableauLobattoIIIC(Float64, s)
TableauLobattoIIIC̄(s) = TableauLobattoIIIC̄(Float64, s)
TableauLobattoIIID(s) = TableauLobattoIIID(Float64, s)
TableauLobattoIIIE(s) = TableauLobattoIIIE(Float64, s)
TableauLobattoIIIF(s) = TableauLobattoIIIF(Float64, s)
TableauLobattoIIIG(s) = TableauLobattoIIIG(Float64, s)
