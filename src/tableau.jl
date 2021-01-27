
using DelimitedFiles
using Markdown
using PrettyTables

import LinearAlgebra: istril


@doc raw"""
Holds the tableau of a Runge-Kutta method.

Parameters:
 * `T`: datatype of coefficient arrays

Fields:
 * `name`: symbolic name of the tableau
 * `o`: order of the method
 * `s`: number of stages
 * `a`: coefficients $a_{ij}$ with $ 1 \le i,j \le s$
 * `b`: weights $b_{i}$  with $ 1 \le i \le s$
 * `c`: nodes $c_{i}$  with $ 1 \le i \le s$
"""
struct Tableau{T}
    name::Symbol
    o::Int
    s::Int

    a::Matrix{T}
    b::Vector{T}
    c::Vector{T}

    â::Matrix{T}
    b̂::Vector{T}
    ĉ::Vector{T}

    function Tableau{T}(name,o,s,a,b,c) where {T}
        @assert s > 0 "Number of stages must be > 0"
        @assert s == size(a,1) == size(a,2) == length(b) == length(c)

        ã = convert(Matrix{T}, a)
        b̃ = convert(Vector{T}, b)
        c̃ = convert(Vector{T}, c)

        â = similar(ã)
        b̂ = similar(b̃)
        ĉ = similar(c̃)

        â .= a .- ã
        b̂ .= b .- b̃
        ĉ .= c .- c̃

        new(name,o,s,ã,b̃,c̃,â,b̂,ĉ)
    end

    function Tableau{T}(name,o,a,b,c) where {T}
        Tableau{T}(name,o,length(c),a,b,c)
    end
end

Tableau(name::Symbol, o::Int, s::Int, a::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}) where {T} = Tableau{T}(name,o,s,a,b,c)
Tableau(name::Symbol, o::Int, a::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}) where {T} = Tableau{T}(name,o,a,b,c)


Base.hash(tab::Tableau, h::UInt) = hash(tab.o, hash(tab.s, hash(tab.a, hash(tab.b, hash(tab.c, hash(tab.â, hash(tab.b̂, hash(tab.ĉ, hash(:Tableau, h)))))))))

Base.:(==)(tab1::Tableau, tab2::Tableau) = (tab1.o == tab2.o
                                         && tab1.s == tab2.s
                                         && tab1.a == tab2.a
                                         && tab1.b == tab2.b
                                         && tab1.c == tab2.c
                                         && tab1.â == tab2.â
                                         && tab1.b̂ == tab2.b̂
                                         && tab1.ĉ == tab2.ĉ)

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

"Read Runge-Kutta tableau from file."
function from_file(dir::AbstractString, name::AbstractString)
    file = string(dir, "/", name, ".tsv")

    # Reads and parse Tableau metadata from file
    f = open(file, "r")
    header = readline(f)
    close(f)

    if header[1] == '#'
        header = split(header[2:end])
    else
        header = ()
    end

    if length(header) ≥ 1
        o = Base.parse(Int, header[1])
    else
        o = 0
    end

    if length(header) ≥ 2
        s = Base.parse(Int, header[2])
    else
        s = 0
    end

    if length(header) ≥ 3
        T = Core.eval(Main, Meta.parse(header[3]))
    else
        T = Float64
    end


    # TODO Read data in original format (e.g., Rational).
    #      For this we need to save tableaus as jld or hdf5.
    # tab_array = readdlm(file, T)
    tab_array = readdlm(file, comments=true)

    if s == 0
        s = size(tab_array, 1) - 1
    end

    @assert s == size(tab_array, 1) - 1 == size(tab_array, 2) - 1

    @info("Reading Runge-Kutta tableau $(name) with $(s) stages and order $(o) from file\n$(file)")

    from_array(tab_array, Symbol(name), o)
end

"Write Runge-Kutta tableau to file."
function to_file(dir::AbstractString, tab::Tableau{T}) where {T}
    # tab_array = zeros(T, S+1, S+1)
    # tab_array[1:S, 2:S+1] = tab.a
    # tab_array[S+1, 2:S+1] = tab.b
    # tab_array[1:S, 1] = tab.c
    # tab_array[S+1, 1] = tab.order

    tab_array = to_array(tab)
    header = string("# ", tab.o, " ", tab.s, " ", T, "\n")
    file   = string(dir, "/", tab.name, ".tsv")

    @info("Writing Runge-Kutta tableau $(tab.name) with $(tab.s) stages and order $(tab.o) to file\n$(file)")

    f = open(file, "w")
    write(f, header)
    writedlm(f, float(tab_array))
    close(f)

    # TODO Write data in original format (e.g., Rational).
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
