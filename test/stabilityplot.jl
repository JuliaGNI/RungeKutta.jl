using CairoMakie
using RungeKutta
using Test

import TOML

const RB = Rational{BigInt}

implicit_euler(::Type{T}) where {T} = Tableau{T}(:implicit_euler, 1, [1;;], [1], [1])

# The children of a stability plot, found by their type.
child(p, ::Type{P}) where {P} = only(filter(c -> c isa P, p.plots))
children(p, ::Type{P}) where {P} = filter(c -> c isa P, p.plots)

# The arguments of a child, as Makie stores them after its conversion.
arguments(c) = (c[1][], c[2][], c[3][])

# Each entry of `a` equals the entry of `b` converted to the element type of `a`, up to 8 ULPs.
function matches(a::AbstractArray{F}, b) where {F}
    bF = F.(b)
    size(a) == size(bF) &&
        all(isequal(u, v) || isapprox(u, v; rtol = 8eps(F)) for (u, v) in zip(a, bF))
end

# The vertices of the boundary line of a stability plot that are not NaN.
function boundary_vertices(p)
    lines = child(child(p, Makie.Contour), Makie.Lines)
    filter(v -> !any(isnan, v), lines[1][])
end

stableplot(tab; kwargs...) = stabilityplot!(Axis(Figure()[1, 1]), tab; kwargs...)

@testset "Weak dependency on Makie" begin
    project = TOML.parsefile(joinpath(pkgdir(RungeKutta), "Project.toml"))
    @test haskey(project["weakdeps"], "Makie")
    @test !haskey(project["deps"], "Makie")
    @test project["extensions"]["RungeKuttaMakieExt"] == "Makie"
    @test haskey(project["compat"], "Makie")
    @test Base.get_extension(RungeKutta, :RungeKuttaMakieExt) isa Module
end

@testset "Stubs and docstrings in the main package" begin
    for s in (:stabilityplot, :stabilityplot!)
        @test s in names(RungeKutta)
        @test haskey(Base.Docs.meta(RungeKutta), Base.Docs.Binding(RungeKutta, s))
    end
end

@testset "Return values and saved images" begin
    fap = stabilityplot(TableauRK4())
    @test fap isa Makie.FigureAxisPlot
    p = stabilityplot!(Axis(Figure()[1, 1]), TableauRK4())
    @test p isa Makie.Plot{RungeKutta.stabilityplot}
    dir = mktempdir()
    for (name, fig) in (("fap.png", fap), ("plot.png", Makie.parent_scene(p)))
        file = joinpath(dir, name)
        save(file, fig)
        @test read(file)[1:4] == [0x89, 0x50, 0x4e, 0x47]
    end
end

@testset "Plotted data equals the region data in $T" for T in (Float64, BigFloat, RB)
    tab = implicit_euler(T)
    kw = (; xlims = (-1.0, 3.0), ylims = (-2.0, 2.0), n = 5)
    r = stability_region(tab; kw...)
    @test eltype(r.absR) == float(T)
    p = stableplot(tab; kw...)
    for c in (child(p, Makie.Contourf), child(p, Makie.Contour))
        x, y, absR = arguments(c)
        @test x isa Vector{Float64}
        @test y isa Vector{Float64}
        @test absR isa Matrix{Float32}
        @test matches(x, r.x)
        @test matches(y, r.y)
        @test matches(absR, r.absR)
        @test isinf.(absR) == isinf.(r.absR)
        @test (x[3], y[3]) == (1.0, 0.0)
        @test isinf(absR[3, 3])
        @test !any(isnan, absR)
    end
end

@testset "The defaults equal those of stability_region" begin
    r = stability_region(TableauRK4())
    x, y, absR = arguments(child(stableplot(TableauRK4()), Makie.Contourf))
    @test matches(x, r.x)
    @test matches(y, r.y)
    @test matches(absR, r.absR)
end

@testset "The keywords reach stability_region" begin
    p = stableplot(TableauRK4(); xlims = (-2.0, 2.0), ylims = (-1.0, 1.0), n = 7)
    x, y, absR = arguments(child(p, Makie.Contourf))
    @test matches(x, collect(range(-2.0, 2.0; length = 7)))
    @test matches(y, collect(range(-1.0, 1.0; length = 7)))
    @test size(absR) == (7, 7)
end

@testset "A change of an attribute recomputes the children" begin
    p = stableplot(TableauRK4(); n = 5)
    @test length(child(p, Makie.Contourf)[1][]) == 5
    p.n[] = 9
    @test length(child(p, Makie.Contourf)[1][]) == 9
    @test size(child(p, Makie.Contour)[3][]) == (9, 9)
end

@testset "The fill is the set where |R| ≤ 1" begin
    dir = mktempdir()
    tab = implicit_euler(Float64)
    for (k, (xlims, ylims, npolygons)) in enumerate((((-6.0, -4.0), (-1.0, 1.0), 1),
        ((1.2, 1.5), (-0.1, 0.1), 0)))
        fig = Figure()
        p = stabilityplot!(Axis(fig[1, 1]), tab; xlims, ylims, n = 21)
        polygons = child(child(p, Makie.Contourf), Makie.Poly)[1][]
        @test length(polygons) == npolygons
        @test isempty(boundary_vertices(p))
        file = joinpath(dir, "window$k.png")
        save(file, fig)
        @test read(file)[1:4] == [0x89, 0x50, 0x4e, 0x47]
    end
end

@testset "The boundary is the set where |R| = 1" begin
    h = 0.04
    vertices = boundary_vertices(stableplot(implicit_euler(Float64)))
    @test length(vertices) ≥ 50
    @test all(abs(hypot(v[1] - 1, v[2]) - 1) ≤ h for v in vertices)
end

@testset "The coordinate axes" begin
    fig = Figure()
    ax = Axis(fig[1, 1])
    p = stabilityplot!(ax, TableauRK4(); xlims = (1.0, 3.0), ylims = (0.5, 2.0))
    hl = children(p, Makie.HLines)
    vl = children(p, Makie.VLines)
    @test length(hl) == 1
    @test length(vl) == 1
    @test only(hl)[1][] == 0
    @test only(vl)[1][] == 0
    reset_limits!(ax)
    limits = ax.finallimits[]
    @test all(>(0), limits.origin)
    @test all(limits.origin .<= (1.0, 0.5))
    @test all(limits.origin .+ limits.widths .>= (3.0, 2.0))
end

@testset "A refused window reaches the caller" begin
    tab = TableauRK4()
    for kw in ((; n = 1), (; xlims = (2.0, 1.0)), (; ylims = (-Inf, 0.0)))
        @test_throws r"^ArgumentError: the window must have finite, increasing limits and n ≥ 2$" stabilityplot(
            tab;
            kw...)
        @test_throws r"^ArgumentError: the window must have finite, increasing limits and n ≥ 2$" stabilityplot!(
            Axis(Figure()[1, 1]), tab; kw...)
    end
end

@testset "A partitioned tableau has no stability plot" begin
    pair = PartitionedTableau(TableauRK4())
    @test_throws r"^MethodError: no method matching stability_region\(" stabilityplot(pair)
end
