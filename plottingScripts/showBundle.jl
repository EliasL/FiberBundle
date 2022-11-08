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

function show_fb(b::FB)
    show_array(b.status, b.L)
end

function show_array(a::AbstractArray, L=nothing)
    if L===nothing
        L = round(Int, sqrt(length(a)))
    end
    m = reshape(a, (L, L))
    show_matrix(m)
end

function show_matrix(m::AbstractMatrix, use_shift=false)
    if use_shift
        m = shift_spanning_cluster(m)
    end
    h = heatmap(m, c=:glasbey_category10_n256, legend=:none, showaxis = false, ticks=false)
    plot(h)
end
