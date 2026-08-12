"""
Compute a vector spanning the one-dimensional nullspace of the square matrix `A`.

Uses a column-pivoted QR factorisation, which unlike the singular value
decomposition is implemented generically in `LinearAlgebra` and thus works for
arbitrary floating point types such as `BigFloat`. With `A * P = Q * R` and `A`
of rank `n-1`, the last diagonal entry of `R` vanishes, so `z = [y; 1]` with
`R[1:n-1,1:n-1] * y = -R[1:n-1,n]` spans the nullspace of `R` and `P * z` that
of `A`. The result is normalised to unit length and its first nonzero entry is
made positive, so that the vector is determined by `A` alone and does not depend
on the pivot order.

`A` must be at least `2×2` and of rank exactly `n-1`, otherwise an `ArgumentError`
is thrown rather than a vector returned that spans nothing. Pivoting orders the
diagonal of `R` by decreasing magnitude, so that is the case precisely when the
last entry has collapsed relative to the first and the one before it has not. For
the matrices this is used on the two are separated by some seventy orders of
magnitude, so the exact threshold does not matter.
"""
function _nullvector(A::AbstractMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(A)

    if n < 2
        throw(ArgumentError("A one-dimensional nullspace needs a matrix of size at least 2×2."))
    end

    F = LinearAlgebra.qr(A, LinearAlgebra.ColumnNorm())
    R = F.R

    tol = n * sqrt(eps(real(float(T)))) * abs(R[1,1])

    if iszero(R[1,1]) || abs(R[n,n]) > tol || abs(R[n-1,n-1]) ≤ tol
        throw(ArgumentError("Matrix is not of rank n-1, so its nullspace is not one-dimensional."))
    end

    y = LinearAlgebra.UpperTriangular(R[1:n-1, 1:n-1]) \ (-R[1:n-1, n])

    w = Vector{T}(undef, n)
    w[F.p] = vcat(y, one(T))
    w = LinearAlgebra.normalize(w)

    return w .* sign(w[findfirst(!iszero, w)])
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
