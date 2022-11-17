using Plots
using DataStructures

include("../burningMan.jl")

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

function shift_spanning_cluster(m::AbstractMatrix)
    shift = get_ideal_shift(m)
    return circshift(m, shift)
end

function shift_spanning_cluster!(b::FB)
    m = reshape(b.status, (b.L, b.L))
    shift = get_ideal_shift(m)
    m = circshift(m, shift)
    b.status = reshape(m, (b.N))
end

function plot_fb(b::FB; show=true, axes=false)
    plot_array(b.status, L=b.L, show=show, axes=axes)
end

function plot_array(a::AbstractArray; L=nothing, show=true, axes=false)
    if L===nothing
        L = round(Int, sqrt(length(a)))
    end
    m = reshape(a, (L, L))
    plot_matrix(m, show=show, axes=axes)
end

function plot_matrix(m::AbstractMatrix; use_shift=false, show=true, axes=false)
    if use_shift
        m = shift_spanning_cluster(m)
    end
    h = heatmap(m, c=:glasbey_category10_n256, legend=:none, aspect_ratio=:equal, showaxis = axes, ticks=axes)
    p = plot(h)
    if show
        display(p)
    end
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

function plot_fb_cm(b::FB)
    return plot!(b.cluster_cm_x, b.cluster_cm_y, seriestype = :scatter)
end