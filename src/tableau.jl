
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


"Print Runge-Kutta coefficients."
function Base.show(io::IO, tab::Tableau)
    print(io, "Runge-Kutta Tableau ", tab.name, " with ", tab.s, " stages and order ", tab.o)
    print(io, "  a = ", tab.a)
    print(io, "  b = ", tab.b)
    print(io, "  c = ", tab.c)
end
