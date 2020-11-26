
import Polynomials
import Polynomials: Polynomial


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

    D(k) = Polynomials.derivative(Polynomial(BigFloat[0, 1])^(k-1) * Polynomial(BigFloat[-1, 1])^k, k-1)
    c = sort(real.(Polynomials.roots(D(s))))
    c[end] = 1; c
end

@doc raw"""
The Radau weights are implicitly given by the so-called simplifying assumption $B(s)$:
```math
\sum \limits_{j=1}^{s} b_{j} c_{j}^{k-1} = \frac{1}{k}  \qquad k = 1 , \, ... , \, s .
```
"""
function get_radau_weights(s)
    if s == 1
        throw(ErrorException("Radau weights for one stage are not defined."))
    end
    solve_simplifying_assumption_b(get_radau_nodes(s))
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
    solve_simplifying_assumption_c(get_radau_nodes(s))
end


"Radau IIA tableau with s stages"
function TableauRadauIIA(s, T=Float64)
    Tableau{T}(Symbol("RadauIIA($s)"), 2s-1, get_radau_coefficients(s), get_radau_weights(s), get_radau_nodes(s))
end
