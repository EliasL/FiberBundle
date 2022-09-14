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

function show_bundle(grid::Vector)
    p = make_grid(grid)
    display(p)
end

function show_bundle(grid1::Vector, grid2::Vector)
    display(
       [ make_grid(grid1),
        make_grid(grid2)]
    ) 
end