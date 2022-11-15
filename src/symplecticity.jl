
function check_symplecticity(tab::Tableau{T}; atol=16eps(T), rtol=16eps(T)) where {T}
    a, b = tab.a, tab.b
    [isapprox(b[i] * a[i,j] + b[j] * a[j,i], b[i] * b[j]; atol=atol, rtol=rtol) for i in axes(a,1), j in axes(a,2)]
end

function issymplectic(tab::Tableau; kwargs...)
    all(check_symplecticity(tab; kwargs...))
end


function check_symplecticity(tab::PartitionedTableau{T}; atol=16eps(T), rtol=16eps(T)) where {T}
    a, b = tab.q.a, tab.q.b
    ā, b̄ = tab.p.a, tab.p.b
    ([isapprox(b[i] * ā[i,j] + b̄[j] * a[j,i], b[i] * b̄[j]; atol=atol, rtol=rtol) for i in axes(a,1), j in axes(a,2)],
     [isapprox(b[i], b̄[i]; atol=atol, rtol=rtol) for i in eachindex(b,b̄)])
end

function issymplectic(tab::PartitionedTableau; kwargs...)
    all(all.(check_symplecticity(tab; kwargs...)))
end


function compute_symplecticity_error(tab::Tableau)
    a, b = tab.a, tab.b
    [b[i] * a[i,j] + b[j] * a[j,i] - b[i] * b[j] for i in axes(a,1), j in axes(a,2)]
end


function get_symplectic_conjugate_coefficients!(a̅::AbstractMatrix{T}, a::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    @assert size(a) == size(a̅)
    @assert length(b) == size(a,1) == size(a,2)

    for i in axes(a̅, 1)
        for j in axes(a̅, 2)
            a̅[i,j] = b[j] / b[i] * ( b[i] - a[j,i] )
        end
    end

    return a̅
end

get_symplectic_conjugate_coefficients(a, b) = get_symplectic_conjugate_coefficients!(zero(a), a, b)
get_symplectic_conjugate_coefficients(tab::Tableau) = Tableau(tab.name, tab.s, get_symplectic_conjugate_coefficients(tab.a, tab.b), tab.b, tab.c)


function symplecticize(tab::Tableau; name=nothing, T=Float64, R∞=tab.R∞)
    a̅ = get_symplectic_conjugate_coefficients(tab.a, tab.b)
    Tableau{T}(name === nothing ? Symbol(tab.name, "S") : name, tab.o, (tab.a .+ a̅) ./ 2, tab.b, tab.c; R∞=R∞)
end
