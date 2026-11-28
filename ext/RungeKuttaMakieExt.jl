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
    x = lift(r -> r.x, region)
    y = lift(r -> r.y, region)
    absR = lift(r -> r.absR, region)
    contourf!(p, x, y, absR; levels = [0.0, 1.0])
    contour!(p, x, y, absR; levels = [1.0])
    hlines!(p, 0; yautolimits = false)
    vlines!(p, 0; xautolimits = false)
    return p
end

end
