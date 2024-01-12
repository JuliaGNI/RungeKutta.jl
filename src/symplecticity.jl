
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


function symplecticity_error(a::AbstractMatrix, b::AbstractVector)
    T = promote_type(eltype(a), eltype(b))
    T[b[i] * a[i,j] + b[j] * a[j,i] - b[i] * b[j] for i in axes(a,1), j in axes(a,2)]
end

symplecticity_error(tab::Tableau) = symplecticity_error(tab.a, tab.b)


function symplectic_conjugate_coefficients(a::AbstractMatrix, b::AbstractVector)
    @assert length(b) == size(a,1) == size(a,2)
    T = promote_type(eltype(a), eltype(b))
    T[b[j] / b[i] * ( b[i] - a[j,i] ) for i in axes(a, 2), j in axes(a, 1)]
end


"""
    SymplecticTableau(tab::Tableau)

Generates a new tableau with symplectizied coefficients.
"""
function SymplecticTableau(tab::Tableau{TT}; name=Symbol("Symplectic", tab.name), T=TT, R∞=tab.R∞) where {TT}
    a = (tab.a .+ symplectic_conjugate_coefficients(tab.a, tab.b)) ./ 2
    Tableau{T}(name, tab.o, a, tab.b, tab.c; R∞=R∞)
end

"""
    SymplecticConjugateTableau(tab::Tableau)

Generates a new tableau with symplectic conjugate coefficients.
"""
function SymplecticConjugateTableau(tab::Tableau{TT}; name=tab.name, T=TT, R∞=tab.R∞) where {TT}
    Tableau{T}(name, tab.s, symplectic_conjugate_coefficients(tab.a, tab.b), tab.b, tab.c; R∞=R∞)
end

"""
    SymplecticPartitionedTableau(tab::Tableau)

Generates a partitioned tableau with tab and ist symplectic adjoint.
"""
SymplecticPartitionedTableau(tab::Tableau{TT}; name=Symbol("Symplectic", tab.name), T=TT, R∞=tab.R∞) where {TT} = PartitionedTableau{T}(name, tab, SymplecticConjugateTableau(tab); R∞=R∞)
