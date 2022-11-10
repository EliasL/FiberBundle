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

function plot_fb(b::FB; show=true)
    plot_array(b.status, L=b.L, show=show)
end

function plot_array(a::AbstractArray; L=nothing, show=true)
    if L===nothing
        L = round(Int, sqrt(length(a)))
    end
    m = reshape(a, (L, L))
    plot_matrix(m, show=show)
end

function plot_matrix(m::AbstractMatrix; use_shift=false, show=true)
    if use_shift
        m = shift_spanning_cluster(m)
    end
    h = heatmap(m, c=:glasbey_category10_n256, legend=:none, showaxis = false, ticks=false)
    p = plot(h)
    if show
        display(p)
    end
end

function plot_fb_axes(b::FB, minor_axes::AbstractVector, major_axes::AbstractVector)
    

function plot_fb_cm(b::FB)
    plot!(b.cluster_cm_x, b.cluster_cm_y, seriestype = :scatter)
end