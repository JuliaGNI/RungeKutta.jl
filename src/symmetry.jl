
function check_symmetry(tab::Tableau{T}; atol=16*eps(T), rtol=16*eps(T)) where {T}
    symmetric = falses(tab.s, tab.s)

    for i in axes(symmetric, 1)
        for j in axes(symmetric, 2)
            symmetric[i,j] = isapprox(tab.a[tab.s+1-i, tab.s+1-j] + tab.a[i,j], tab.b[j], atol=atol, rtol=rtol)
        end
    end

    symmetric
end

function issymmetric(tab::Tableau{T}; atol=16*eps(T), rtol=16*eps(T)) where {T}
    all(check_symmetry(tab; atol=atol, rtol=rtol))
end
