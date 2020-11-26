
using LinearAlgebra
using Polynomials
using SpecialPolynomials


function get_gauss_nodes(s)
    p = zeros(BigFloat, s+1); p[end] = 1
    c = sort(real.(roots(convert(Polynomial, ShiftedLegendre(p)))))
end

function get_gauss_weights(s)
    c = get_gauss_nodes(s)
    p = zeros(BigFloat, s+1); p[end] = 1
    P = convert(Polynomial, ShiftedLegendre(p))
    D = Polynomials.derivative(P)
    
    inti(i) = begin
        I = Polynomials.integrate( ( P ÷ Polynomial([-c[i], 1]) )^2 )
        I(1) - I(0)
    end
    
    b = [ inti(i) / D(c[i])^2  for i in 1:s ]
end

function vandermonde_matrix_inverse(x::Vector{T}) where {T}
    local n = length(x)

    local L::Matrix{T} = zeros(T,n,n)
    local U::Matrix{T} = Matrix{T}(I, n, n)
    local V::Matrix{T}

    L[1,1] = 1
    for i in 2:n
        for j in 1:i
            p = 1
            for k in 1:i
                if k ≠ j
                    p *= (x[j] - x[k])
                end
            end
            L[i,j] = 1/p
        end
    end

    i = 1
    for j in i+1:n
        U[i,j] = - U[i,j-1] * x[j-1]
    end

    for i in 2:n
        for j in i+1:n
            U[i,j] = U[i-1,j-1] - U[i,j-1] * x[j-1]
        end
    end

    V = *(U,L)
end


function TableauGauss(s::Int, T=Float64)

    # order
    o = 2s

    # obtain Gauss-Legendre nodes and weights
    b = get_gauss_weights(s)
    c = get_gauss_nodes(s)

    # create Lagrange polynomial
    vdm = vandermonde_matrix_inverse(c)

    # compute monomial basis functions and corresponding integrals
    poly_ints = []
    for i in 1:s
        y = zeros(s)
        y[i] = 1
        mon = *(vdm, y)
        push!(poly_ints, Polynomials.integrate(Polynomials.Polynomial(mon)))
    end

    # compute Runge-Kutta coefficients
    a = zeros(s,s)
    for i in 1:s
        for j in 1:s
            a[i,j] = poly_ints[j](c[i])
        end
    end

    Tableau{T}(Symbol("Gauss", s), o, a, b, c)
end
