using RungeKutta
using Documenter
using Weave


module RungeKuttaWeaves
    using Markdown
    using Markdown: MD, Paragraph, LineBreak
    using Plots
    using PrettyTables
    using RungeKutta
    using RungeKutta: to_array
    using RungeKutta: get_gauss_nodes, get_gauss_weights, get_gauss_coefficients
    using RungeKutta: get_lobatto_nodes, get_lobatto_weights,
                      get_lobatto_coefficients_a, get_lobatto_coefficients_b,
                      get_lobatto_coefficients_c, get_lobatto_coefficients_c̄,
                      get_lobatto_coefficients_f
    using RungeKutta: get_radau_1_nodes, get_radau_1_weights, get_radau_1_coefficients,
                      get_radau_2_nodes, get_radau_2_weights, get_radau_2_coefficients
    using SymPy
    using SymPy: latex

    "Markdown-print Runge-Kutta tableau with SymPy coefficients."
    function Base.show(io::IO, ::MIME"text/markdown", tab::Tableau{Sym})
        show(io, "text/markdown", Markdown.parse("Runge-Kutta Tableau $(tab.name) with $(tab.s) stages and order $(tab.o):"))

        tab_arr = simplify.(to_array(tab))
        str_arr = latex.(tab_arr)
        str_arr[tab.s+1,1] = ""

        strio = IOBuffer()
        pretty_table(strio, str_arr,
                        backend = :latex,
                        vlines = [1],
                        hlines = [tab.s],
                        noheader = true)
        tab_latex = String(take!(strio))

        tab_markdown = replace(tab_latex, "tabular" => "array")
        tab_markdown = replace(tab_markdown, "\\begin{table}" => "```math")
        tab_markdown = replace(tab_markdown, "\\end{table}" => "```")

        print(io, tab_markdown * "\n")
    end
end


weave("src/gauss.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)

weave("src/radauIA.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)

weave("src/radauIIA.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)

weave("src/lobatto.jmd",
         out_path = "src",
         doctype = "github",
         mod = RungeKuttaWeaves)


makedocs(;
    modules=[RungeKutta],
    authors="Michael Kraus",
    repo="https://github.com/JuliaGNI/RungeKutta.jl/blob/{commit}{path}#L{line}",
    sitename="RungeKutta.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliagni.github.io/RungeKutta.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Runge-Kutta Methods" => "rungekutta.md",
        "Gauß Methods" => "gauss.md",
        "Radau IA Methods" => "radauIA.md",
        "Radau IIA Methods" => "radauIIA.md",
        "Lobatto III Methods" => "lobatto.md",
        "Tableau" => "tableau.md",
        "Diagnostics" => "diagnostics.md",
        "Library" => "library.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGNI/RungeKutta.jl",
)
