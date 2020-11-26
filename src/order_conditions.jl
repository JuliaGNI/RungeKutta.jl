
function check_order_conditions_b(tab::Tableau{T}, k) where {T}
    local res::T = 0

    for i in 1:tab.s
        res += tab.b[i] * tab.c[i]^(k-1)
    end

    return isapprox(res, 1/k, atol=16*eps(T), rtol=16*eps(T))
end


function check_order_conditions_c(tab::Tableau{T}, k) where {T}
    local order  = falses(tab.s)
    local res::T

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
    local res::T

    for i in axes(tab.a, 1)
        res = 0
        for j in axes(tab.a, 2)
            res += tab.b[j] * tab.c[j]^(k-1) * tab.a[j,i]
        end
        order[i] = isapprox(res, tab.b[i] * (1 - tab.c[i]^k)/k, atol=16*eps(T), rtol=16*eps(T))
    end

    order
end
