using RungeKutta
using Documenter
using DocumenterCitations
using Weave

bib = CitationBibliography(joinpath(@__DIR__, "RungeKutta.bib"))

module RungeKuttaWeaves
    using Markdown
    using Markdown: MD, Paragraph, LineBreak
    using Plots
    using PrettyTables
    using RungeKutta
    using RungeKutta.Tableaus: get_gauss_nodes, get_gauss_weights, get_gauss_coefficients
    using RungeKutta.Tableaus: get_lobatto_nodes, get_lobatto_weights,
                               get_lobatto_a_coefficients, get_lobatto_b_coefficients,
                               get_lobatto_c_coefficients, get_lobatto_c̄_coefficients,
                               get_lobatto_d_coefficients, get_lobatto_e_coefficients,
                               get_lobatto_f_coefficients, get_lobatto_g_coefficients
    using RungeKutta.Tableaus: get_radau_1_nodes, get_radau_1_weights, get_radau_1_coefficients,
                               get_radau_2_nodes, get_radau_2_weights, get_radau_2_coefficients

    import SymPy
    import SymPy: latex, simplify

    """
        symtype()

    Return `Sym{T}` for `T` being the underlying type of `Sym(1)`.
    """
    symtype() = typeof(SymPy.Sym(1))

    
    "Markdown-print Runge-Kutta tableau with SymPy coefficients."
    function Base.show(io::IO, ::MIME"text/markdown", tab::Tableau{symtype()})
        show(io, "text/markdown", Markdown.parse("Runge-Kutta Tableau $(tab.name) with $(tab.s) stages and order $(tab.o):"))

        tab_arr = simplify.(convert(Matrix, tab))
        str_arr = latex.(tab_arr)
        str_arr[tab.s+1,1] = ""

        strio = IOBuffer()
        pretty_table(strio, LatexCell.(str_arr),
                        backend = Val(:latex),
                        vlines = [1],
                        hlines = [tab.s],
                        show_header = false)
        tab_latex = String(take!(strio))

        tab_markdown = replace(tab_latex, "tabular" => "array")
        # tab_markdown = replace(tab_markdown, "\\begin{table}" => "```math")
        # tab_markdown = replace(tab_markdown, "\\end{table}" => "```")
        tab_markdown = "```math\n" * tab_markdown * "```\n"

        print(io, tab_markdown * "\n")
    end
end


weave("src/gauss.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)

weave("src/radau1.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)

weave("src/radau2.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)

weave("src/lobatto.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)


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
        "Bibliography" => "bibliography.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaGNI/RungeKutta.jl",
    devurl = "latest",
    devbranch = "main",
)
