
@doc raw"""
The Gauss coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function gauss_coefficients(::Type{T}, s) where {T}
    solve_simplifying_assumption_c(gauss_legendre_nodes(T,s))
end

gauss_coefficients(s) = gauss_coefficients(BigFloat, s)


reference(::Val{:Gauss}) = """
References:

    John C. Butcher.
    Implicit Runge-Kutta processes.
    Mathematics of Computation, Volume 18, Pages 50-64, 1964.
    doi: 10.1090/S0025-5718-1964-0159424-9.

    John C. Butcher.
    Gauss Methods. 
    In: Engquist B. (eds). Encyclopedia of Applied and Computational Mathematics. Springer, Berlin, Heidelberg. 2015.
    doi: 10.1007/978-3-540-70529-1_115.
"""

@doc raw"""
Gauss tableau with s stages

```julia
TableauGauss(::Type{T}, s)
TableauGauss(s) = TableauGauss(Float64, s)
```
The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

The nodes and weights are those of the Gauss-Legendre quadrature rule, i.e. the roots of the
shifted Legendre polynomial $P_s (2x-1)$ and the corresponding interpolatory weights, taken
from `QuadratureRules.gauss_legendre_nodes` and `QuadratureRules.gauss_legendre_weights`. The
coefficients follow from the simplifying assumption $C(s)$, cf. [`gauss_coefficients`](@ref).
Prescribing no node leaves all $2s$ parameters free, which is what gives the method its
order $2s$.

""" * reference(Val(:Gauss)) function TableauGauss(::Type{T}, s) where {T}
    Tableau{T}(:Gauss, 2s, gauss_coefficients(s), gauss_legendre_weights(BigFloat, s),
               gauss_legendre_nodes(BigFloat, s); R∞=(-1)^s)
end

TableauGauss(s) = TableauGauss(Float64, s)
