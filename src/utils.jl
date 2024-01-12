"Legendre polynomial P_s(x) of degree s defined on the interval [-1..+1]."
function _legendre(j::Integer, p::Polynomial{T}) where {T}
    if j ≤ 0
        return Polynomial(one(T))
    elseif j == 1
        return p
    else
        return ( (2j-1) * p * _legendre(j-1, p) - (j-1) * _legendre(j-2, p) ) / j
    end
end

"Legendre polynomial of degree s shifted to the interval [0..1], i.e., P_s(2x-1)."
function _shifted_legendre(s, T=BigFloat)
    _legendre(s, Polynomial(T[-1, 2]))
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
