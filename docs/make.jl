using RungeKutta
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "RungeKutta.bib"))

"Helpers for the `@eval` blocks of the pages with symbolic tableaus."
module SymbolicTableaus
using Markdown
using PrettyTables
using RungeKutta

import SymPyPythonCall
import SymPyPythonCall: latex, simplify

"""
    symtype()

Return `Sym{T}` for `T` being the underlying type of `Sym(1)`.
"""
symtype() = typeof(SymPyPythonCall.Sym(1))

"Markdown-print Runge-Kutta tableau with SymPy coefficients."
function Base.show(io::IO, ::MIME"text/markdown", tab::Tableau{symtype()})
    show(io,
        "text/markdown",
        Markdown.parse("Runge-Kutta Tableau $(tab.name) with $(tab.s) stages and order $(tab.o):"))

    tab_arr = simplify.(convert(Matrix, tab))
    str_arr = latex.(tab_arr)
    str_arr[tab.s + 1, 1] = ""

    strio = IOBuffer()
    pretty_table(strio, LatexCell.(str_arr),
        backend = :latex,
        table_format = RungeKutta.butcher_latex_tableau_format(tab),
        show_column_labels = false,
        show_row_number_column = false
    )
    tab_latex = String(take!(strio))

    tab_markdown = replace(tab_latex, "tabular" => "array")
    # tab_markdown = replace(tab_markdown, "\\begin{table}" => "```math")
    # tab_markdown = replace(tab_markdown, "\\end{table}" => "```")
    tab_markdown = "```math\n" * tab_markdown * "```\n"

    print(io, tab_markdown * "\n")
end

"Markdown of the tableaus `f(s)` for `s` in `stages`, as the value of an `@eval` block."
function tableaus_markdown(f, stages)
    Markdown.parse(join(repr("text/markdown", f(s)) for s in stages))
end
end

makedocs(
    sitename = "RungeKutta.jl",
    authors = "Michael Kraus",
    plugins = [bib],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules = [RungeKutta],
    pages = [
        "Home" => "index.md",
        "Runge-Kutta Methods" => "rungekutta.md",
        "Tableau" => "tableau.md",
        "Diagnostics" => "diagnostics.md",
        "Tabulated Methods" => "methods.md",
        "Gauß Methods" => "gauss.md",
        "Radau I Methods" => "radau1.md",
        "Radau II Methods" => "radau2.md",
        "Lobatto III Methods" => "lobatto.md",
        "Library" => "library.md",
        "Bibliography" => "bibliography.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaGNI/RungeKutta.jl",
    devurl = "latest",
    devbranch = "main"
)
