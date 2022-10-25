using Luxor
using JLD2
using MathTeXEngine
using LaTeXStrings
using Plots
using Random

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



function draw_seeds(L, color_stress)

    global_path = "data/"
    ts = search_for_t(global_path)
    NRS = ["UNR", "SNR", "CNR"]
    distribution(t, NR) = "t=$t Uniform $NR"

    seed = 9
    ps = 10 #Pixel size

    f(t, NR) = load(full_name(global_path, L, distribution(t, NR)))

    t_settings = [0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.7, 0.8, 0.9]
    layout = (length(NRS),length(t_settings))
    lx = layout[1] #layout x
    ly = layout[2] #layout y
    key = color_stress ? "spanning_cluster_tension" : "spanning_cluster_state"
    #println(f(0.0, NR))
    grids = reshape([reshape(f(t, NR)["$key/$seed"], (L, L)) for NR=NRS, t=t_settings], (lx,ly))
    kn(t, NR) = round(f(t, NR)["spanning_cluster_step/$seed"]/(L*L), digits=2)
    grid_names = reshape([latexstring("$NR, \$t=$t\$") for NR=NRS, t=t_settings], (lx,ly))
    pixel_L = L*ps
    title_space = ceil(Int, pixel_L/2)
    subtitle_space = ceil(Int, pixel_L/5)
    spacing_x = ceil(Int, pixel_L/20)
    spacing_y = spacing_x + subtitle_space
    image_size_x = pixel_L*lx+(lx-1)*spacing_x
    image_size_y = pixel_L*ly+(ly)*spacing_y + title_space
    font_size = 4*pixel_L/32
    title_font_size = font_size * 4/3

    stress =  color_stress ? " stress" : ""
    Drawing(image_size_x, image_size_y, "plots/Visualizations/Spanning/Spanning clusters L=$L$stress.png")
    background(1,1,1) #White background
    fontface("Computer Modern")
    fontsize(title_font_size)
    Luxor.text(latexstring("Spanning clusters, \$ L=$L\$"), image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
    fontsize(font_size)
    for j=1:ly, i=1:lx
        origin((pixel_L+spacing_x)*(i-1), (pixel_L+spacing_y)*(j-1) + title_space)
        drawmatrix(grids[i,j], color_stress, ps)
        Luxor.text(grid_names[i, j], pixel_L/2, pixel_L+subtitle_space, halign=:center, valign=:bottom)
    end

    finish()
end

draw_seeds(32, true)
draw_seeds(32, false)
draw_seeds(64, true)
draw_seeds(64, false)