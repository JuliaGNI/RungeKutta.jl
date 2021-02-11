
@doc raw"""
The Gauss nodes are given by the roots of the shifted Legendre polynomial
$P_s (2x-1)$ with $s$ the number of stages.
"""
function get_gauss_nodes(::Type{T}, s) where {T}
    sort(T.(Polynomials.roots(_shifted_legendre(s,T))))
end

get_gauss_nodes(s) = get_gauss_nodes(BigFloat, s)


@doc raw"""
The Gauss weights are given by the following integrals
```math
b_i = \bigg( \frac{dP}{dx} (c_i) \bigg)^{-2} \int \limits_0^1 \bigg( \frac{P(x)}{x - c_i} \bigg)^2 dx ,
```
where $P(x)$ denotes the shifted Legendre polynomial
$P(x) = P_s (2x-1)$ with $s$ the number of stages.
"""
function get_gauss_weights(::Type{T}, s) where {T}
    c = get_gauss_nodes(T,s)
    P = _shifted_legendre(s,T)
    D = Polynomials.derivative(P)
    
    inti(i) = begin
        I = Polynomials.integrate( ( P ÷ Polynomial(T[-c[i], 1]) )^2 )
        I(1) - I(0)
    end
    
    b = [ inti(i) / D(c[i])^2  for i in 1:s ]
end

get_gauss_weights(s) = get_gauss_weights(BigFloat, s)


@doc raw"""
The Gauss coefficients are implicitly given by the so-called simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function get_gauss_coefficients(::Type{T}, s) where {T}
    solve_simplifying_assumption_c(get_gauss_nodes(T,s))
end

get_gauss_coefficients(s) = get_gauss_coefficients(BigFloat, s)


"""
Gauss tableau with s stages

```julia
TableauGauss(::Type{T}, s)
TableauGauss(s) = TableauGauss(Float64, s)
```

The constructor takes the number of stages `s` and optionally the element type `T` of the tableau.

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
function TableauGauss(::Type{T}, s) where {T}
    Tableau{T}(Symbol("Gauss($s)"), 2s, get_gauss_coefficients(s), get_gauss_weights(s), get_gauss_nodes(s))
end

TableauGauss(s) = TableauGauss(Float64, s)
