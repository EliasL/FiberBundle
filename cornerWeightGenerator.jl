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

function reshape_neighbours(l)
    m = reshape(l, (3,3))
    m[3,3] = m[2,3]
    m[2,3] = m[1,3]
    m[1,3] = m[3,2]
    m[3,2] = m[2,2]
    m[2,2] = 2
    return m
end

function generate_neighbours()
    nr_neighbours = 257
    neighbours = [zeros(Int, (3, 3)) for i in 1:nr_neighbours]
    for i in 1:nr_neighbours
        n = reverse(parse.(Int, split(string(i-1, base=2, pad=9), "")))
        neighbours[i] = reshape_neighbours(n)
    end
    return neighbours
end


function remove_symmetries(n)
    for grid in n
        syms = [grid, rotr90(grid), rot180(grid), rotl90(grid), reverse(grid, dims=1), reverse(grid, dims=2), rotr90(reverse(grid, dims=1)), rotr90(reverse(grid, dims=2))]
        duplicates = findall(x->x in syms,n)[2:end]
        deleteat!(n, duplicates)
    end
    return n
end

function compute_strength(m::Matrix)
    # Here are the strength values
    # 121
    # 2_2
    # 121
    strength = sum([m[2,1], m[1,2], m[3,2], m[2,3]])*2
    strength += sum([m[1,1], m[3,1], m[1,3], m[3,3]])

    if m == rot180(m)
        strength += 0.5
    end
    return strength
end

neighbours = generate_neighbours()

# Remove symmetries
neighbours = remove_symmetries(neighbours)

# Sort by strength 
neighbours = sort(neighbours, by=compute_strength)

extra_strength = [
    0,0,0,0,0, # 5
    0,0,0,0,0, # 10
    0,0,0,0,0, # 15
    0,0,0,0,0, # 20
    0,0,0,0,0, # 25
    0,0,0,0,0, # 30
    0,0,0,0,0, # 35
    0,0,0,0,0, # 40
    0,0,0,0,0, # 45
    0,0,0,0,0, # 50
    0,0        # 52
]

n = length(neighbours)

L = 3
columns = 5
layout = (columns,ceil(Int64, n/columns))
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
text("Neighbours", image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
fontsize(font_size)

for j=1:ly, i=1:lx
    neighbour_index = i+(j-1)*lx
    if neighbour_index<=length(neighbours)
        origin((L+spacing_x)*(i-1), (L+spacing_y)*(j-1) + title_space)
        drawmatrix(neighbours[neighbour_index])
        id = i+(j-1)*lx
        s = compute_strength(neighbours[neighbour_index])
        ss = extra_strength[neighbour_index]
        text("$id: $s + $ss", L/2, L+subtitle_space, halign=:center, valign=:bottom)
    end
end

finish()