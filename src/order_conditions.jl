
function check_order_conditions_b(coeff::Tableau{T}, k) where {T}
    local res::T = 0

    for i in 1:coeff.s
        res += coeff.b[i] * coeff.c[i]^(k-1)
    end

    return isapprox(res, 1/k, atol=16*eps(T), rtol=16*eps(T))
end


function check_order_conditions_c(coeff::Tableau{T}, k) where {T}
    local order  = falses(coeff.s)
    local res::T

    for i in 1:size(coeff.a, 1)
        res = 0
        for j in 1:size(coeff.a, 2)
            res += coeff.a[i,j] * coeff.c[j]^(k-1)
        end
        order[i] = isapprox(res, coeff.c[i]^k/k, atol=16*eps(T), rtol=16*eps(T))
    end

    order
end


function check_order_conditions_d(coeff::Tableau{T}, k) where {T}
    local order  = falses(coeff.s)
    local res::T

    for i in 1:size(coeff.a, 1)
        res = 0
        for j in 1:size(coeff.a, 2)
            res += coeff.b[j] * coeff.c[j]^(k-1) * coeff.a[j,i]
        end
        order[i] = isapprox(res, coeff.b[i] * (1 - coeff.c[i]^k)/k, atol=16*eps(T), rtol=16*eps(T))
    end

    order
end
