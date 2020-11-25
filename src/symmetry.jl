
function check_symmetry(coeff::Tableau{T}) where {T}
    symmetric = falses(coeff.s, coeff.s)

    for i in 1:size(coeff.a, 1)
        for j in 1:size(coeff.a, 2)
            symmetric[i,j] = isapprox(coeff.a[coeff.s+1-i, coeff.s+1-j] + coeff.a[i,j], coeff.b[j], atol=16*eps(T), rtol=16*eps(T))
        end
    end

    symmetric
end
