"""
Compute a vector spanning the one-dimensional nullspace of the square matrix `A`.

Uses a column-pivoted QR factorisation, which unlike the singular value
decomposition is implemented generically in `LinearAlgebra` and thus works for
arbitrary floating point types such as `BigFloat`. With `A * P = Q * R` and `A`
of rank `n-1`, the last diagonal entry of `R` vanishes, so `z = [y; 1]` with
`R[1:n-1,1:n-1] * y = -R[1:n-1,n]` spans the nullspace of `R` and `P * z` that
of `A`. The result is normalised to unit length.
"""
function _nullvector(A::AbstractMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(A)
    F = LinearAlgebra.qr(A, LinearAlgebra.ColumnNorm())
    R = F.R

    y = LinearAlgebra.UpperTriangular(R[1:n-1, 1:n-1]) \ (-R[1:n-1, n])

    w = Vector{T}(undef, n)
    w[F.p] = vcat(y, one(T))

    return LinearAlgebra.normalize(w)
end


function istriustrict(A::AbstractMatrix)
    m, n = size(A)
    if m == n == 1
        if A[1,1] ≠ 0
            return false
        end
    else
        @inbounds for j in 1:min(n,m-1), i in j:m
            if A[i,j] ≠ 0
                return false
            end
        end
    end
    return true
end

function istrilstrict(A::AbstractMatrix)
    m, n = size(A)
    if m == n == 1
        if A[1,1] ≠ 0
            return false
        end
    else
        @inbounds for j in 2:n, i in 1:min(j,m)
            if A[i,j] ≠ 0
                return false
            end
        end
    end
    return true
end
