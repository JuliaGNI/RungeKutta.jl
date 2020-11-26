
function check_order_conditions_b(tab::Tableau{T}, k) where {T}
    local res::T = 0

    for i in 1:tab.s
        res += tab.b[i] * tab.c[i]^(k-1)
    end

    return isapprox(res, 1/k, atol=16*eps(T), rtol=16*eps(T))
end


function check_order_conditions_c(tab::Tableau{T}, k) where {T}
    local order  = falses(tab.s)
    local res::T = 0

    for i in axes(tab.a, 1)
        res = 0
        for j in axes(tab.a, 2)
            res += tab.a[i,j] * tab.c[j]^(k-1)
        end
        order[i] = isapprox(res, tab.c[i]^k/k, atol=16*eps(T), rtol=16*eps(T))
    end

    order
end


function check_order_conditions_d(tab::Tableau{T}, k) where {T}
    local order  = falses(tab.s)
    local res::T = 0

    for i in axes(tab.a, 1)
        res = 0
        for j in axes(tab.a, 2)
            res += tab.b[j] * tab.c[j]^(k-1) * tab.a[j,i]
        end
        order[i] = isapprox(res, tab.b[i] * (1 - tab.c[i]^k)/k, atol=16*eps(T), rtol=16*eps(T))
    end

    order
end


@doc raw"""
Compute the weights by solving the simplifying assumption $B(s)$:
```math
\sum \limits_{j=1}^{s} b_{j} c_{j}^{k-1} = \frac{1}{k}  \qquad k = 1 , \, ... , \, s .
```
"""
function solve_simplifying_assumption_b(c)
    s = length(c)
    M = [ c[j]^(k-1) for k in 1:s, j in 1:s ]
    r = [ 1 / k for k in big.(1:s) ]
    M \ r
end

@doc raw"""
Compute the coefficients by solving the simplifying assumption $C(s)$:
```math
\sum \limits_{j=1}^{s} a_{ij} c_{j}^{k-1} = \frac{c_i^k}{k}  \qquad i = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function solve_simplifying_assumption_c(c)
    s = length(c)
    M = [ c[j]^(k-1) for k in 1:s, j in 1:s ]
    
    row(i) = begin
        r = [ c[i]^k / k for k in 1:s ]
        M \ r
    end
    
    vcat([row(i)' for i in 1:s]...)
end

@doc raw"""
Compute the coefficients by solving the simplifying assumption $D(s)$:
```math
\sum \limits_{i=1}^{s} b_i c_{i}^{k-1} a_{ij} = \frac{b_j}{k} ( 1 - c_j^k)  \qquad j = 1 , \, ... , \, s , \; k = 1 , \, ... , \, s .
```
"""
function solve_simplifying_assumption_d(b,c)
    @assert axes(b) == axes(c)

    s = length(c)
    M = [ b[i] * c[i]^(k-1) for k in 1:s, i in 1:s ]
    
    row(j) = begin
        r = [ b[j] / k * (1 - c[j]^k) for k in 1:s ]
        M \ r
    end
    
    hcat([row(j) for j in 1:s]...)
end
