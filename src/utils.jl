
_big(x) = x
_big(x::Number) = big(x)
_big(x::String) = parse(BigFloat, x)

function _big(x::Expr)
    y = x
    y.args .= _big.(y.args)
    return  y
end

macro big(x)
    return esc(_big(x))
end


function istriustrict(A::AbstractMatrix)
    m, n = size(A)
    if m == n == 1
        if A[1,1] != 0
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
        if A[1,1] != 0
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
