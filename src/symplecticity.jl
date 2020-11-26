
function check_symplecticity(tab::Tableau{T}; atol=16*eps(T), rtol=16*eps(T)) where {T}
    symplectic = falses(tab.s, tab.s)

    for i in axes(tab.a, 1)
        for j in axes(tab.a, 2)
            symplectic[i,j] = isapprox(tab.b[i] * tab.a[i,j] + tab.b[j] * tab.a[j,i], tab.b[i] * tab.b[j], atol=atol, rtol=rtol)
        end
    end

    symplectic
end

function issymplectic(tab::Tableau{T}; atol=16*eps(T), rtol=16*eps(T)) where {T}
    all(check_symplecticity(tab; atol=atol, rtol=rtol))
end


function compute_symplecticity_error(tab::Tableau)
    [tab.b[i] * tab.a[i,j] + tab.b[j] * tab.a[j,i] - tab.b[i] * tab.b[j] for i in axes(tab.a, 1), j in axes(tab.a,2)]
end


function get_symplectic_conjugate_coefficients!(a::AbstractMatrix{T}, b::AbstractVector{T}, a̅::AbstractMatrix{T}) where {T}
    @assert size(a) == size(a̅)
    @assert length(b) == size(a,1) == size(a,2)

    for i in axes(a̅, 1)
        for j in axes(a̅, 2)
            a̅[i,j] = b[j] / b[i] * ( b[i] - a[j,i] )
        end
    end

    return a̅
end

get_symplectic_conjugate_coefficients(a, b) = get_symplectic_conjugate_coefficients!(a, b, zero(a))
get_symplectic_conjugate_coefficients(tab::Tableau) = Tableau(tab.name, tab.s, get_symplectic_conjugate_coefficients(tab.a, tab.b), tab.b, tab.c)


function symplecticize(tab::Tableau; name=nothing, T=Float64)
    a̅ = get_symplectic_conjugate_coefficients(tab.a, tab.b)
    Tableau{T}(name === nothing ? Symbol(tab.name, "S") : name, tab.o, (tab.a .+ a̅) ./ 2, tab.b, tab.c)
end

