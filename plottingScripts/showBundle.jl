using Plots
using DataStructures

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

function shift_spanning_cluster!(b::FB)
    # Get shift
    m = reshape(b.status, (b.L, b.L))
    shift = get_ideal_shift(m)

    # Shift status
    m = circshift(m, shift)
    b.status = reshape(m, (b.N))

    # Shift cm
    b.cluster_cm_x = mod1.(b.cluster_cm_x.+shift[2], b.L)
    b.cluster_cm_y = mod1.(b.cluster_cm_y.+shift[1], b.L)
end

function plot_fb(b::FB; show=true, axes=false, use_shift=true)
    L=b.L
    spanning=b.spanning_cluster_id

    if use_shift
        shift_spanning_cluster!(b)
    end
    m = reshape(b.status, (L, L))
    L = size(m,1)
    nr_clusters = maximum(m)
    nr_colors = nr_clusters
    background_color = :black
    spanning_color = :red
    colors = vcat([background_color],[RGBA(rand(3)...) for _ in eachindex(1:nr_colors)])
    if spanning != -1
        colors[spanning+1] = spanning_color
    end

    c = cgrad(colors)

    # We only want positive states, ie, clusters and 0 for
    # all others
    clamp!(m, 0, Inf)
    image_size = maximum([500,L])+100
    h = heatmap(m, c=c, legend=:none, aspect_ratio=:equal, bg_inside = nothing,
    showaxis = axes, ticks=axes, size=(image_size, image_size))
    p = plot(h)
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
    for (c, r) in R
        cmx = b.cluster_cm_x[c]
        cmy = b.cluster_cm_y[c]
        plot!(circleShape(cmx, cmy, r), seriestype=[:shape,],lw=0.5, c=:blue,
        linecolor=:black, legend=false, fillalpha = 0.2, aspect_ratio=1)
    end
end

function plot_fb_cm(b::FB)
    return plot!(b.cluster_cm_x, b.cluster_cm_y, seriestype = :scatter)
end

function generate_illustrations()
    
    path = "data/"
    t = 0.0
    L=64
    α = 2.0
    seed = 1
    nr = "LLS"
    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed")
    nr="CLS"
    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed")
    L=1024
    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed")
    nr="LLS"
    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed")
end

function save_picture(L, nr, t, α, seed, name, path="data/")
    settings = make_settings("Uniform", L, t, nr, α, path)
    b = get_bundle_from_settings(settings, seed=seed)
    p = plot_fb(b, show=false)
    # We always save the plot as latest_plot, so we can just copy that file
    cp("latest_plot.png", "plots/Visualizations/differenceIllustrations/$name.png", force=true)

end

function test(seeds=1)
    nr = "CLS"
    path = "data/"
    t = 0.1
    L=32
    α = 2.0
    settings = make_settings("Uniform", L, t, nr, α, path)
    bundles = get_bundles_from_settings(settings, seeds=seeds, step=-0)
    for b in bundles
        p = plot_fb(b, show=false)
        display(p)
    end
end

#generate_illustrations()
#test(1:10)