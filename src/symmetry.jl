
function check_symmetry(tab::Tableau{T}; atol=16eps(T), rtol=16eps(T)) where {T}
    a, b, s = tab.a, tab.b, tab.s
    [isapprox(a[s+1-i, s+1-j] + a[i,j], b[j]; atol=atol, rtol=rtol) for i in axes(a,1), j in axes(a,2)]
end

function issymmetric(tab::Tableau{T}; kwargs...) where {T}
    all(check_symmetry(tab; kwargs...))
end


function check_symmetry(tab::PartitionedTableau{T}; kwargs...) where {T}
    (check_symmetry(tab.q; kwargs...),
     check_symmetry(tab.p; kwargs...))
end

function issymmetric(tab::PartitionedTableau{T}; kwargs...) where {T}
    issymmetric(tab.q) && issymmetric(tab.p)
end
