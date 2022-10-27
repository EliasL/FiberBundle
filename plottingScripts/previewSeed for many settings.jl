using Luxor
using JLD2
using MathTeXEngine
using LaTeXStrings
using Plots
using Random
using DataStructures

include("../support/dataManager.jl")

function drawmatrix(A::Matrix, color_stress=true, pixel_size = 1)
    L = size(A,1)
    tiles = Tiler(L*pixel_size, L*pixel_size, L, L, margin=0)
    cur_colors = palette(:glasbey_category10_n256)
    nr_colors = 256
    if color_stress
        A = ceil.(Int, A./ maximum(A) .* (nr_colors-1))
    end
    stress_colors = cgrad(:heat, nr_colors, categorical=true)
    for (pos, n) in tiles
        pos = pos .+ L*pixel_size/2
        fiber = A[tiles.currentcol, tiles.currentrow]
        if color_stress
            box(pos, tiles.tilewidth, tiles.tileheight, :clip)
            if fiber == 0
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


function get_ideal_shift(m)
    # Shift_matrix tries to shift the matrix so that a cluster does#
    # not cross the periodic boarder. In other words, as few broken fibers along
    # the boarder as possible
    L = size(m,1)

    #Find the largest cluster 
    count = counter(filter(f -> f>0, m))
    largest_cluster = collect(keys(count))[argmax(collect(values(count)))]
    part_of_largest_cluster = reshape( m .== largest_cluster, (L,L))

    # Now we check which row and which column has the fewest broken fibers
    # row and col might be swapped here. Haven't chekced.
    min_row = argmin([sum(view(part_of_largest_cluster, i, :)) for i in 1:L])
    min_col = argmin([sum(view(part_of_largest_cluster, :, i)) for i in 1:L])
    return (L-min_row, L-min_col)
end


function load_file(α, NR, settings)
    if NR=="UNR"
        α=0.0
    end
    settings = filter(s -> s["a"]==α && s["nr"]==NR, settings)
    @assert length(settings) < 2 "There are multiple possibilities"
    @assert length(settings) != 0 "There is no file maching these settings α=$α nr=$NR"
    setting = settings[1]
    return load(get_file_name(setting))
end

function draw_seeds(L)

    path = "data/"
    dist = "Uniform"
    settings = search_for_settings(path, dist)
    # We now have all the settings, but we only want to use some of them
    settings =  filter(
        s -> s["t"] == 0 &&
             s["L"] == L,
             settings
    )

    seed = 9
    ps = 10 #Pixel size

    # This function specifies an α and nr and reutrns the file with this setting
    NRS = ["UNR", "SNR", "CNR"]
    α_settings = [1.0, 1.5, 2.0, 3.0, 5.0, 9.0, 15.0]
    lx, ly = (length(NRS),length(α_settings))
    # This will store the best shift of the matrix to center the largest cluster
    ideal_shifts = []


    # we first draw the cluster image since we need those to calculate the ideal shifts
    for key in ["spanning_cluster_state", "spanning_cluster_tension"]
        grids = reshape([reshape(load_file(α, NR, settings)["$key/$seed"], (L, L)) for NR=NRS, α=α_settings], (lx,ly))
        grid_names = reshape([latexstring("$NR, \$a=$α\$") for NR=NRS, α=α_settings], (lx,ly))
        pixel_L = L*ps
        title_space = ceil(Int, pixel_L/2)
        subtitle_space = ceil(Int, pixel_L/5)
        spacing_x = ceil(Int, pixel_L/20)
        spacing_y = spacing_x + subtitle_space
        image_size_x = pixel_L*lx+(lx-1)*spacing_x
        image_size_y = pixel_L*ly+(ly)*spacing_y + title_space
        font_size = 4*pixel_L/32
        title_font_size = font_size * 4/3

        
        stress =  key=="spanning_cluster_tension" 
        Drawing(image_size_x, image_size_y, "plots/Visualizations/Spanning/Spanning clusters L=$L"*"$(stress ? " stress" : "").png")
        background(1,1,1) #White background
        fontface("Computer Modern")
        fontsize(title_font_size)
        Luxor.text(latexstring("Spanning clusters, \$ L=$L\$"), image_size_x/2, title_space*2/3, halign=:center, valign=:bottom)
        fontsize(font_size)
        for j=1:ly, i=1:lx
            origin((pixel_L+spacing_x)*(i-1), (pixel_L+spacing_y)*(j-1) + title_space)

            # Get the shift of the matrix and store it so we can use it when drawing stress
            m = grids[i,j]
            if key=="spanning_cluster_state"
                shift = get_ideal_shift(m)
                push!(ideal_shifts, shift)
            else
                shift = popfirst!(ideal_shifts)
            end

            drawmatrix(circshift(m, shift), stress, ps)
            Luxor.text(grid_names[i, j], pixel_L/2, pixel_L+subtitle_space, halign=:center, valign=:bottom)
        end

        finish()
    end


end

draw_seeds(32)
draw_seeds(64)