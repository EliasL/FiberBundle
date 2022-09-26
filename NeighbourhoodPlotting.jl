using Luxor
using JLD2
using MathTeXEngine
using LaTeXStrings

include("NeighbourhoodWeighting.jl")

function drawmatrix(A::Matrix)
    L = size(A,1)
    tiles = Tiler(L, L, L, L, margin=0)
    
    for (pos, n) in tiles
        pos = pos .+ L/2
        fiber = A[ceil(Int, pos[1]), ceil(Int, pos[2])]
        if fiber == 1
            box(pos, tiles.tilewidth*0.8, tiles.tileheight*0.8, :clip)
            if mod(tiles.currentcol+tiles.currentrow, 2) == 1
                background(0.2,0.2,0.3)
            else
                background(0.2,0.2,0.1)
            end
        elseif fiber == 2
            box(pos, tiles.tilewidth*0.8, tiles.tileheight*0.8, :clip)
            background(0.8,0.2,0.2)
        else 
            box(pos, tiles.tilewidth*0.8, tiles.tileheight*0.8, :clip)
            background(0.9,0.9,0.9)
        end
        clipreset()
    end
end


allNeighbours = generate_neighbours()

# Remove symmetries
neighbours = remove_symmetries(allNeighbours)

# Sort by strength 
neighbours = sort(neighbours, by=compute_strength)
allNeighboursSorted = sort(allNeighbours, by=compute_strength)

function plot_neighbours(neighbours, title, path, columns=5)

    n = length(neighbours)

    L = 3
    layout = (columns,ceil(Int64, n/columns))
    lx = layout[1] #layout x
    ly = layout[2] #layout y

    title_space = L/2
    subtitle_space = L/5
    spacing_x = L/3
    spacing_y = spacing_x + subtitle_space
    image_size_x = L*lx+(lx-1)*spacing_x
    image_size_y = L*ly+(ly)*spacing_y + title_space
    font_size = 4*L/32
    title_font_size = font_size * 4/3

    Drawing(image_size_x, image_size_y, path*title*".pdf")
    fontface("Computer Modern")
    fontsize(title_font_size)
    text(title, image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
    fontsize(font_size)

    for j=1:ly, i=1:lx
        neighbour_index = i+(j-1)*lx
        # 
        if neighbour_index<=length(neighbours)
            origin((L+spacing_x)*(i-1), (L+spacing_y)*(j-1) + title_space)
            drawmatrix(neighbours[neighbour_index])
            id = i+(j-1)*lx
            s = compute_strength(neighbours[neighbour_index])
            text("$id: $s", L/2, L+subtitle_space, halign=:center, valign=:bottom)
        end
    end
    finish()
end

plot_neighbours(neighbours, "Neighbours", "plots/")
plot_neighbours(allNeighbours, "All Neighbours", "plots/")
plot_neighbours(allNeighboursSorted, "All Neighbours Sorted", "plots/")