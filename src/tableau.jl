
using LinearAlgebra
using Markdown
using PrettyTables

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
eachstage(tab::Tableau) = 1:tab.s
coefficients(tab::Tableau) = tab.a
weights(tab::Tableau) = tab.b
nodes(tab::Tableau) = tab.c

isexplicit(tab::Tableau) = istrilstrict(tab.a)
isimplicit(tab::Tableau) = !isexplicit(tab)
isdiagnonallyimplicit(tab::Tableau) = tab.s != 1 && !istrilstrict(tab.a) && istril(tab.a)
isfullyimplicit(tab::Tableau) = (tab.s == 1 && tab.a[1,1] != 0) || (!istrilstrict(tab.a) && !istril(tab.a))


const tf_butcher_tableau = TextFormat(
    up_right_corner     = ' ',
    up_left_corner      = ' ',
    bottom_left_corner  = ' ',
    bottom_right_corner = ' ',
    up_intersection     = ' ',
    left_intersection   = ' ',
    right_intersection  = ' ',
    middle_intersection = ' ',
    bottom_intersection = ' ',
    column              = '│',
    row                 = ' '
)

"Pretty-print Runge-Kutta tableau."
function Base.show(io::IO, tab::Tableau)
    print(io, "\nRunge-Kutta Tableau $(tab.name) with $(tab.s) stages and order $(tab.o):\n")
    arr = convert(Array{Any}, to_array(tab))
    arr[tab.s+1,1] = ""
    pretty_table(io, arr,
                    tf = tf_butcher_tableau,
                    vlines = [1],
                    body_hlines = [tab.s],
                    body_hlines_format = ('─','┼','─','─'),
                    equal_columns_width = true,
                    noheader = true)
end

"Markdown-print Runge-Kutta tableau."
function Base.show(io::IO, ::MIME"text/markdown", tab::Tableau)
    show(io, "text/markdown", Markdown.parse("Runge-Kutta Tableau $(tab.name) with $(tab.s) stages and order $(tab.o):"))

    tab_arr = convert(Array{Any}, to_array(tab))
    tab_arr[tab.s+1,1] = ""

    strio = IOBuffer()
    pretty_table(strio, tab_arr,
                    backend = :latex,
                    vlines = [1],
                    hlines = [tab.s],
                    noheader = true)
    tab_latex = String(take!(strio))

    tab_markdown = replace(tab_latex, "tabular" => "array")
    tab_markdown = replace(tab_markdown, "\\begin{table}" => "```math")
    tab_markdown = replace(tab_markdown, "\\end{table}" => "```")

    print(io, tab_markdown)
end
