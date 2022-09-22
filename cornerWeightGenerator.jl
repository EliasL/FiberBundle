using Luxor
using JLD2
using MathTeXEngine
using LaTeXStrings

function drawmatrix(A::Matrix)
    L = size(A,1)
    tiles = Tiler(L, L, L, L, margin=0)
    
    for (pos, n) in tiles
        pos = pos .+ L/2
        fiber = A[ceil(Int, pos[1]), ceil(Int, pos[2])]
        if fiber == 1
            box(pos, tiles.tilewidth*0.8, tiles.tileheight*0.8, :clip)
            background(0.2,0.2,0.2)
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

function reshape_neighbours(l)
    m = reshape(l, (3,3))
    m[1,1] = m[2,2]
    m[2,2] = 2
    return m
end

function generate_neighbours()
    nr_neighbours = 256
    neighbours = [zeros(Int, (3, 3)) for i in 1:256]
    for i in 1:nr_neighbours
        n = parse.(Int, split(string(i, base=2, pad=9), ""))
        neighbours[i] = reshape_neighbours(n)
    end
    return neighbours
end


function remove_rotations(n)
    for grid in n
        rotations = [rotr90(grid), rot180(grid), rotl90(grid)]
        for rotation in rotations
            duplicates = findall
            deleteat!(n, findall(x->x==rotation,n))
        end
        push!(n, grid)
    end

    return n
end

function remove_mirror(n)
    for grid in n
        mirrors = [reverse(grid, dims=1), reverse(grid, dims=2)]
        for mirror in mirrors
            deleteat!(n, findall(x->x==mirror,n))
        end
        push!(n, grid)
    end
    return n
end

neighbours = generate_neighbours()

# Remove symmetries
neighbours = remove_rotations(neighbours)
println(length(neighbours))
neighbours = remove_mirror(neighbours)
println(length(neighbours))
# After removing symetries, there are 39 left, so we want to 
# create a even number to display in two colunms

push!(neighbours,zeros(Int64, (3,3)))

L = 3
layout = (2,20)
lx = layout[1] #layout x
ly = layout[2] #layout y
area = lx*ly

title_space = L/2
subtitle_space = L/5
spacing_x = L/3
spacing_y = spacing_x + subtitle_space
image_size_x = L*lx+(lx-1)*spacing_x
image_size_y = L*ly+(ly)*spacing_y + title_space
font_size = 4*L/32
title_font_size = font_size * 4/3

Drawing(image_size_x, image_size_y, "plots/Neighbours.pdf")
fontface("Computer Modern")
fontsize(title_font_size)
text(latexstring("Neighbours"), image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
fontsize(font_size)

for j=1:ly, i=1:lx
    origin((L+spacing_x)*(i-1), (L+spacing_y)*(j-1) + title_space)
    drawmatrix(neighbours[i+(j-1)*lx])
    id = i+(j-1)*lx
    text("$id", L/2, L+subtitle_space, halign=:center, valign=:bottom)
end

finish()