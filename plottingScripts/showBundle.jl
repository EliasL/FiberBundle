using Plots
#using DataStructures

include("../burningMan.jl")
include("../support/dataManager.jl")

function circleShape(x, y, r)
    θ = LinRange(0, 2*π, 100)
    return x .+ r*sin.(θ), y .+ r*cos.(θ)
end

function get_ideal_shift(m::AbstractMatrix)
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

function shift_to_cm(b::FB)
    cmx = b.cluster_cm_x[b.spanning_cluster_id]
    cmy = b.cluster_cm_y[b.spanning_cluster_id]
    return ((round(Int, b.L/2-cmy)), (round(Int, b.L/2-cmx)))
end

function shift_spanning_cluster!(b::FB, cm_shift=true)
    # Get shift
    m = reshape(b.status, (b.L, b.L))
    if cm_shift
        shift = shift_to_cm(b)
    else
        shift = get_ideal_shift(m)
    end
    # Shift status
    m = circshift(m, shift)
    b.status = reshape(m, (b.N))
    # Shift cm
    b.cluster_cm_x = mod1.(b.cluster_cm_x.+shift[2], b.L)
    b.cluster_cm_y = mod1.(b.cluster_cm_y.+shift[1], b.L)
end


function plot_fb(b::FB; show=true, axes=false, use_shift=true, stress=false, cm_shift=true)
    L=b.L
    spanning=argmax(b.cluster_size)
    
    if use_shift
        shift_spanning_cluster!(b, cm_shift)
    end
    if stress
        m = reshape(b.tension, (L,L))
        c =:thermal
    else
        m = reshape(b.status, (L, L))
        clamp!(m, 0, Inf)
        nr_clusters = maximum(m)
        nr_colors = nr_clusters
        background_color = RGBA(0,0,0)
        spanning_color = RGBA(0.9,0.1,0.1)
        colors = vcat([background_color],[RGBA((rand(3))...) for i in 1:nr_colors])
        if spanning != -1
            colors[spanning+1] = spanning_color
        end
        c = colors
    end

    # We only want positive states, ie, clusters and 0 for
    # all others
    image_size = maximum([500,L])+100

    img = map(x -> c[x+1], m)
    p = plot(img, legend=:none, aspect_ratio=:equal, bg_inside = nothing,
    showaxis = axes, ticks=axes, size=(image_size, image_size))

    if show
        display(p)
    end
    save("latest_plot.png", p)
    return p
end

function plot_fb_axes(b::FB, minor_axes::AbstractMatrix, major_axes::AbstractMatrix,
                             minor_values::AbstractVector, major_values::AbstractVector)

    lines_x = []
    lines_y = []
    maj_w = sqrt.(abs.(major_values)) / maximum(sqrt.(abs.(major_values))) * b.L/4
    min_w = sqrt.(abs.(minor_values)) / maximum(sqrt.(abs.(major_values))) * b.L/4
    for c in 1:b.c
        x,y = b.cluster_cm_x[c], b.cluster_cm_y[c]
        for slope in [minor_axes[c, :]*min_w[c], major_axes[c, :]*maj_w[c]]
            points_x = []
            points_y = []
            p1 = (x+slope[1], y+slope[2])
            p2 = (x-slope[1], y-slope[2])
            push!(points_x, p1[1])
            push!(points_x, p2[1])
            push!(points_y, p1[2])
            push!(points_y, p2[2])
            push!(lines_x, points_x)
            push!(lines_y, points_y)
            width = log(b.cluster_size[c]) / log(maximum(b.cluster_size)) * 4
            plot!(points_x, points_y, width=width)
        end
    end
    #return plot!(lines_x, lines_y, width=hcat(maj_width, maj_width))
end

function plot_gyration_radi(b::FB, R; nr=:all)
    if nr != :all
        R = sort(collect(enumerate(R)),by= x -> x[2], rev=true)[1:nr]
    else
        R = enumerate(R)
    end
    s = b.cluster_size[b.spanning_cluster_id]
    r_ = round(R[1][2]; digits=2)
    for (c, r) in R
        cmx = b.cluster_cm_x[c]
        cmy = b.cluster_cm_y[c]
        plot!(circleShape(cmx, cmy, r), seriestype=[:shape,],lw=0.5, c=:blue,
        linecolor=:black, legend=false, fillalpha = 0.2, aspect_ratio=1, title="r:$r_, s:$s")
    end
end

function plot_fb_cm(b::FB)
    return plot!(b.cluster_cm_x, b.cluster_cm_y, seriestype = :scatter)
end

function save_picture(L, nr, t, α, seed, name, save_path="", data_path="data/")
    settings = make_settings(L, t, nr, α, data_path)
    b = get_bundles_from_settings(settings, seeds=seed, spanning=true)
    p = plot_fb(b, show=false)
    # We always save the plot as latest_plot, so we can just copy that file
    if save_path != ""
        cp("latest_plot.png", "$save_path/$name.png", force=true)
    end
end

function test(seeds=1)
    nr = "CLS"
    t = 0.1
    L=32
    α = 2.0
    settings = make_settings(L, t, nr, α)
    bundles = get_bundles_from_settings(settings, seeds=seeds, step=-0)
    for b in bundles
        p = plot_fb(b, show=false)
        display(p)
    end
end

#test(1:10)