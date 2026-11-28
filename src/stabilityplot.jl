"""
    stabilityplot(tab::Tableau; xlims = (-6.0, 2.0), ylims = (-4.0, 4.0), n = 201)

Plot the stability region of `tab` in the complex plane: the set where `|R(z)| ≤ 1` is filled,
its boundary `|R(z)| = 1` is drawn as a line, and the real and imaginary axes are drawn as two
lines through 0. The data come from [`stability_region`](@ref), with the same keywords and
defaults, so a plot costs `n²` evaluations of the stability function. A window or a grid that
`stability_region` refuses raises its `ArgumentError`.

`stabilityplot` is a Makie recipe and has a method only when Makie is loaded. A Makie backend,
such as CairoMakie or GLMakie, loads Makie and draws the plot. `stabilityplot` returns a
`Makie.FigureAxisPlot`, and takes the usual `figure` and `axis` keywords of Makie:

```julia
using CairoMakie
stabilityplot(TableauGauss(2); axis = (; aspect = DataAspect(), xlabel = "Re z", ylabel = "Im z"))
```

See [`stabilityplot!`](@ref) to add the plot to an existing axis.
"""
function stabilityplot end

"""
    stabilityplot!(ax, tab::Tableau; xlims = (-6.0, 2.0), ylims = (-4.0, 4.0), n = 201)

Add the stability region of `tab` to the axis `ax`, as [`stabilityplot`](@ref) draws it, and
return the plot object. A change of `xlims`, `ylims` or `n` on the plot object recomputes the
region. A Makie backend must be loaded to draw the plot.
"""
function stabilityplot! end
