using Luxor
using JLD2
using MathTeXEngine
using LaTeXStrings

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*".jld2"
file() = load(full_name(global_path, L, distribution))

function drawmatrix(A::Matrix)
    L = size(A,1)
    tiles = Tiler(L, L, L, L, margin=0)
    for (pos, n) in tiles
        pos = pos .+ L/2
        fiber = A[ceil(Int, pos[1]), ceil(Int, pos[2])]
        if fiber == -1
            ngon(pos, tiles.tilewidth/1.6, 6, 0, :clip)
            background(0.2,0.2,0.2)
        elseif fiber == -3
            box(pos, tiles.tilewidth*0.6, tiles.tileheight*0.6, :clip)
            background(0.2,0.2,0.2)
        else 
            #box(pos, tiles.tilewidth*0.8, tiles.tileheight*0.8, :clip)
            #background(0.15,0.15,0.15)
        end
        clipreset()
    end
end

global_path = "data/"
distribution = "Uniform"
L = 128
seed = 1
f = file()
layout = (5,1)
lx = layout[1] #layout x
ly = layout[2] #layout y
area = lx*ly
nr_stored_states = 9 # Check the save file
states = [floor(Int,nr_stored_states/area*i) for i in 1:area][1:area]

grids = reshape([reshape(f["sample_states/$seed"][i, :], (L, L)) for i in states], (lx,ly))
title_space = L/2
subtitle_space = L/5
spacing_x = L/20
spacing_y = spacing_x + subtitle_space
image_size_x = L*lx+(lx-1)*spacing_x
image_size_y = L*ly+(ly)*spacing_y + title_space
font_size = 4*L/32
title_font_size = font_size * 4/3

Drawing(image_size_x, image_size_y, "plots/$distribution sample view.pdf")
fontface("Computer Modern")
fontsize(title_font_size)
text(latexstring("$distribution distribution, \$ L=$L\$"), image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
fontsize(font_size)
for j=1:ly, i=1:lx
    origin((L+spacing_x)*(i-1), (L+spacing_y)*(j-1) + title_space)
    drawmatrix(grids[i,j])
    kN = states[i+(j-1)*lx]
    text(latexstring("k/N = 0.$kN"), L/2, L+subtitle_space, halign=:center, valign=:bottom)
end

finish()