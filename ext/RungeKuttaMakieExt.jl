module RungeKuttaMakieExt

using Makie: Makie, @recipe, contourf!, contour!, hlines!, vlines!, lift
using RungeKutta: stability_region

import RungeKutta: stabilityplot, stabilityplot!

# The defaults of the attributes are those of `stability_region`.
@recipe StabilityPlot (tab,) begin
    "The real limits of the window."
    xlims = (-6.0, 2.0)
    "The imaginary limits of the window."
    ylims = (-4.0, 4.0)
    "The number of grid points in each direction."
    n = 201
end

function Makie.plot!(p::StabilityPlot)
    region = lift(p[1], p.xlims, p.ylims, p.n) do tab, xlims, ylims, n
        stability_region(tab; xlims, ylims, n)
    end
    # The three arguments of the children change in one update of the compute graph; with one
    # observable each, a child can resolve with the new `x` and the old `absR`. `add_input!` is
    # ComputePipeline's, reached through Makie, which does not export it.
    Makie.add_input!(p.attributes, :region, region)
    map!(r -> (r.x, r.y, r.absR), p.attributes, :region, [:x, :y, :absR])
    contourf!(p, p.x, p.y, p.absR; levels = [0.0, 1.0])
    # A fixed colour keeps the line away from the colormap, whose range comes from the data and is
    # NaN on an all-NaN grid.
    contour!(p, p.x, p.y, p.absR; levels = [1.0], color = :black)
    hlines!(p, 0; yautolimits = false)
    vlines!(p, 0; xautolimits = false)
    return p
end

# An `Axis` reads the autolimit flags only of its top-level plots, and the default bounding box of
# a recipe joins all of its children. This bounding box leaves out the children that set either
# flag to false, the coordinate axes, so that a window that excludes 0 keeps its own limits.
function Makie.boundingbox(p::StabilityPlot, space::Symbol = :data)
    bb = Ref(Makie.Rect3d())
    for c in p.plots
        if Makie.to_value(get(c, :xautolimits, true)) &&
           Makie.to_value(get(c, :yautolimits, true))
            Makie.update_boundingbox!(bb, Makie.boundingbox(c, space))
        end
    end
    return bb[]
end

end
