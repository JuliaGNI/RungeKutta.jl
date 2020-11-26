
"Holds the tableau of a Runge-Kutta method."
struct Tableau{T}
    name::Symbol
    o::Int
    s::Int
    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    function Tableau{T}(name,o,s,a,b,c) where {T}
        @assert s > 0 "Number of stages must be > 0"
        @assert s == size(a,1) == size(a,2) == length(b) == length(c)
        new(name,o,s,a,b,c)
    end

    function Tableau{T}(name,o,a,b,c) where {T}
        Tableau{T}(name,o,length(c),a,b,c)
    end
end

Tableau(name::Symbol, o::Int, s::Int, a::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}) where {T} = Tableau{T}(name,o,s,a,b,c)
Tableau(name::Symbol, o::Int, a::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}) where {T} = Tableau{T}(name,o,length(c),a,b,c)


Base.hash(tab::Tableau, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.a, hash(tab.b, hash(tab.c, hash(:Tableau, h))))))

Base.:(==)(tab1::Tableau, tab2::Tableau) = (tab1.o == tab2.o
                                         && tab1.s == tab2.s
                                         && tab1.a == tab2.a
                                         && tab1.b == tab2.b
                                         && tab1.c == tab2.c)

Base.isapprox(tab1::Tableau, tab2::Tableau; kwargs...) = (
                                            tab1.o == tab2.o
                                         && tab1.s == tab2.s
                                         && isapprox(tab1.a, tab2.a; kwargs...)
                                         && isapprox(tab1.b, tab2.b; kwargs...)
                                         && isapprox(tab1.c, tab2.c; kwargs...))
 
Base.isequal(tab1::Tableau{T1}, tab2::Tableau{T2}) where {T1,T2} = (tab1 == tab2 && T1 == T2 && tab1.name == tab2.name && typeof(tab1) == typeof(tab2))

Base.eltype(::Tableau{T}) where {T} = T

function from_array(arr::AbstractMatrix{T}, name::Symbol, o::Int) where {T}
    @assert size(arr,1) == size(arr,2)

    local s = size(arr,1)-1
    local a = copy(arr[1:s, 2:s+1])
    local b = copy(arr[s+1, 2:s+1])
    local c = copy(arr[1:s, 1    ])

    Tableau{T}(name, o, a, b, c)
end

function to_array(tab::Tableau{T}) where {T}
    local s = tab.s
    local arr = zeros(T, s+1, s+1)
    arr[1:s, 1    ] .= tab.c
    arr[s+1, 2:s+1] .= tab.b
    arr[1:s, 2:s+1] .= tab.a
    return arr
end

Base.convert(::Type{Tableau}, x::AbstractMatrix; name::Symbol=:nonamespecified, o::Int=0) = from_array(x, name, o)
Base.convert(::Type{Matrix{T}}, tab::Tableau) where {T} = convert(Matrix{T}, to_array(tab))

name(tab::Tableau) = tab.name
order(tab::Tableau) = tab.o
nstages(tab::Tableau) = tab.s
coefficients(tab::Tableau) = tab.a
weights(tab::Tableau) = tab.b
nodes(tab::Tableau) = tab.c


"Print Runge-Kutta coefficients."
function Base.show(io::IO, tab::Tableau)
    print(io, "Runge-Kutta Tableau ", tab.name, " with ", tab.s, " stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
end
