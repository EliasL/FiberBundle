using Plots
using Random

plotlyjs()

function make_grid(grid::Vector)
    L = round(Int, sqrt(length(grid)))
    g = reshape(grid, (L,L))
    p = heatmap(g,
            aspect_ratio=:equal,
            colorbar=false)
    return p
end

function show(grid::Vector)
    p = make_grid(grid)
    display(p)
end

function show(grid1::Vector, grid2::Vector)
    display(
       [ make_grid(grid1),
        make_grid(grid2)]
    ) 
end

using Interact, Plots
## Interact.WebIO.install_jupyter_nbextension() # might be helpful if you see `WebIO` warnings in Jupyter
@manipulate throttle=.05 for λ=0:.1:5, μ=0:.1:5
    xs = range(0.0, 1.0, length = 100)
    Plots.plot(xs, x -> λ*x^2 + μ)
end