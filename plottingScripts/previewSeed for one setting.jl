using Luxor
using JLD2
using MathTeXEngine
using LaTeXStrings
using Plots
using Random
using ProgressMeter

include("../support/dataManager.jl")

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*"_bulk.jld2"

function drawmatrix(A::Matrix, color_stress=true, pixel_size = 1)
    L = size(A,1)
    tiles = Tiler(L*pixel_size, L*pixel_size, L, L, margin=0)
    cur_colors = palette(:glasbey_category10_n256)
    nr_colors = 256
    if color_stress
        A = floor.(Int, A./ maximum(A) .* (nr_colors-1)) .+ 1
    end
    stress_colors = cgrad(:heat, nr_colors, categorical=true)
    for (pos, n) in tiles
        pos = pos .+ L*pixel_size/2
        fiber = A[tiles.currentcol, tiles.currentrow]
        if color_stress
            box(pos, tiles.tilewidth, tiles.tileheight, :clip)
            if fiber == 1
                star(pos, tiles.tilewidth*0.5, 4, 0.2, pi/4, action=:clip)
                background(0.2,0.2,0.2)
            else
                background(stress_colors[fiber].r, stress_colors[fiber].g, stress_colors[fiber].b)
            end
        else
                
            if fiber == -1
                box(pos, tiles.tilewidth*0.6, tiles.tileheight*0.6, :clip)
                background(0.2,0.2,0.2)
                
            elseif fiber == -3
                #ngon(pos, tiles.tilewidth/1.6, 6, 0, :clip)
                box(pos, tiles.tilewidth, tiles.tileheight, :clip)
                background(0.2,0.2,0.2)
            else
                box(pos, tiles.tilewidth, tiles.tileheight, :clip)
                background(cur_colors[mod1(fiber, 256)].r, cur_colors[mod1(fiber, 256)].g, cur_colors[mod1(fiber, 256)].b)

            end
        end
        clipreset()
    end
end

function draw_seeds(setting, seed)
    ps = 10 #Pixel size

    color_stress = false
    f = load_file(setting, -1, false)
    layout = (3,3)
    lx = layout[1] #layout x
    ly = layout[2] #layout y
    area = lx*ly
    nr_stored_states = 9 # Check the save file
    states = [floor(Int,nr_stored_states/area*i) for i in 1:area][1:area]
    key = color_stress ? "tension" : "sample_states"
    grids = reshape([reshape(f["$key/$seed"][i, :], (L, L)) for i in states], (lx,ly))
    pixel_L = L*ps
    title_space = ceil(Int, pixel_L/2)
    subtitle_space = ceil(Int, pixel_L/5)
    spacing_x = ceil(Int, pixel_L/20)
    spacing_y = spacing_x + subtitle_space
    image_size_x = pixel_L*lx+(lx-1)*spacing_x
    image_size_y = pixel_L*ly+(ly)*spacing_y + title_space
    font_size = 4*pixel_L/32
    title_font_size = font_size * 4/3

    Drawing(image_size_x, image_size_y, "plots/Visualizations/Progressions/$(setting["name"]) sample view.png")
    background(1,1,1) #White background
    fontface("Computer Modern")
    fontsize(title_font_size)
    Luxor.text(latexstring("$(setting["name"]) distribution, \$ L=$L\$"), image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
    fontsize(font_size)
    for j=1:ly, i=1:lx
        origin((pixel_L+spacing_x)*(i-1), (pixel_L+spacing_y)*(j-1) + title_space)
        drawmatrix(grids[i,j], color_stress, ps)
        kN = states[i+(j-1)*lx]
        Luxor.text(latexstring("k/N = 0.$kN"), pixel_L/2, pixel_L+subtitle_space, halign=:center, valign=:bottom)
    end

    finish()
end

NRS = ["SNR"]
global_path = "data/"
ts = [0.1]
L=512
α = 2.0
seed = 1
p = Progress(length(NRS)*length(ts))
for NR in NRS
    for t in ts
        setting = make_settings("Uniform", L, t, NR, α, global_path)
        draw_seeds(setting, seed)
        ProgressMeter.next!(p)
    end
end
