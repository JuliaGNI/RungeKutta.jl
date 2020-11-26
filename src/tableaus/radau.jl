
import LinearAlgebra
import Polynomials


@doc raw"""
The s-stage Radau nodes are defined as the roots of the following polynomial of degree $s$:
```math
\frac{d^{s-2}}{dx^{s-2}} \big( (x - x^2)^{s-1} \big) .
```
"""
function get_radau_nodes(s)
    if s == 1
        throw(ErrorException("Radau nodes for one stage are not defined."))
    end

    D(k) = Polynomials.derivative(Polynomials.Polynomial(BigFloat[0, 1])^(k-1) * Polynomials.Polynomial(BigFloat[-1, 1])^k, k-1)
    c = sort(real.(Polynomials.roots(D(s))))
    c[end] = 1; c
end

@doc raw"""
The Radau weights can be explicitly computed by the formula
```math
b_j = \frac{1}{s (s-1) P_{s-1}(2 c_j - 1)^2} \qquad j = 1 , \, ... , \, s ,
```
where $P_k$ is the $k$th Legendre polynomial, given by
```math
P_k (x) = \frac{1}{k! 2^k} \big( \frac{d^k}{dx^k} (x^2 - 1)^k \big) .
```
"""
function get_radau_weights(s)
    if s == 1
        throw(ErrorException("Radau weights for one stage are not defined."))
    end

    c = get_radau_nodes(s)
    M = [ c[j]^(k-1) for k in 1:s, j in 1:s ]
    r = [ 1 / k for k in big.(1:s) ]
    M \ r
end

@doc raw"""
The Radau IIA coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_radau_coefficients(s)
    if s == 1
        throw(ErrorException("Radau IIA coefficients for one stage are not defined."))
    end

    c = get_radau_nodes(s)
    M = [ c[j]^(k-1) for k in 1:s, j in 1:s ]
    
    row(i) = begin
        r = [ c[i]^k / k for k in 1:s ]
        M \ r
    end
    
    vcat([row(i)' for i in 1:s]...)
end


"Radau IIA tableau with s stages"
function TableauRadauIIA(s, T=Float64)
    Tableau{T}(Symbol("RadauIIA($s)"), 2s-1, get_radau_coefficients(s), get_radau_weights(s), get_radau_nodes(s))
end
