using RungeKutta
using Documenter

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
        "Library" => "library.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGNI/RungeKutta.jl",
)
