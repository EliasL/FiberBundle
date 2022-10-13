using BenchmarkTools
using Random
using ProgressMeter
using SparseArrays

include("support/neighbourhoodWeighting.jl")

"""
Simulate breakage of an LxL fiber matrix
"""

# Define constants
ALIVE = -1 #::Int64 = -1 # A fiber that has not yet been broken
CURRENT_BORDER = -2 #::Int64 = -2 # The border of the current cluster being explored
PAST_BORDER = -3#::Int64 = -3 # The border of a cluster that has been explored
BROKEN = 0#::Int64 = 0 # A broken fiber that has not been explored
# Status larger than 0 represents the cluster id that the fiber belongs to

 #
#x#
 #
NEIGHBOURS(i, L) = [
    i+L, # Right
    i-L,  # Left
    i-1, # Up
    i+1, # Down
]

###
#x#
###
# In order to optmize the code, we want the neighbourhood to be in 
# this order:
#   147
#   258
#   369
# Note that unlike NEIGHBOURS, NEIGHBOURHOOD includes itself
NEIGHBOURHOOD(i, L) = [
    i-1-L, # Up Left
    i-L,   # Left
    i+1-L, # Down Left
    i-1,   # Up
    i+1,   # Down
    i-1+L, # Up Right
    i+L,   # Right
    i+1+L, # Down Right
]

function fillAdjacent(L::Int64, adjacent_indexes::Function=NEIGHBOURS)
    N = L*L
    nr_adjacent = length(adjacent_indexes(1,L))
    adjacent = zeros(Int64, N, nr_adjacent)
    for i in 1:N
        adjacent[i, :] = findAdjacent(i, L, adjacent_indexes(i, L))
    end
    return adjacent
end

function findAdjacent(i::Int64, L::Int64, adjacent::Vector{Int64})
    # Finds adjacent fibers in periodic boundry contidions
    # The structure of the matrix will be
    # 1 4 6
    # 2   7
    # 3 5 8
    # If the checks are confusing, use the gid above to
    # see what happens.
    # mod1(7,3) = mod1(1,3) = 1
    # mod1(3,3) = mod1(9,3) = 3
    N = L*L
    # Fill in adjacent
    a = adjacent

    if length(a) == 4
        # Right check
        if i>N-L
            a[1] -= N
        end
        # Left check
        if i<=L
            a[2] += N
        end
        # Up check
        if mod1(i, L) == 1
            a[3] += L
        end
        # Down check
        if mod1(i, L) == L
            a[4] -= L
        end
    elseif length(a) == 8
        # Up check
        if mod1(i, L) == 1
            a[1] += L
            a[4] += L
            a[6] += L
        end
        # Down check
        if mod1(i, L) == L
            a[3] -= L
            a[5] -= L
            a[8] -= L
        end
        # Right check
        if i>N-L
            a[6] -= N
            a[7] -= N
            a[8] -= N
        end
        # Left check
        if i<=L
            a[1] += N
            a[2] += N
            a[3] += N
        end
    end

    return a 
end


function resetClusters(status::Vector{Int64}, σ::Vector{Float64})
    # After having explored and assigned numbers to 
    # all the fibers indicating their state, we now 
    # want to reset them to either be BROKEN or ALIVE
    @inbounds @simd for i in eachindex(status)
        s = status[i]

        # The statuses of the fibers will indicate
        # what cluster they belong to. If this fiber
        # has had a status larger than 0, it means that
        # it has belonged to a cluster. 
        if s > 0#BROKEN
            # In this case, we know that it is broken
            # and should be reset
            status[i] = 0#BROKEN

        elseif s < 0#BROKEN
            # This means that this fiber was set as either
            # ALIVE, CURRENT_BORDER or PAST_BORDER.
            # Now we reset it to being just ALIVE
            status[i] = -1 #ALIVE
            # We have already added stress to the boarder, so now
            # AFTER having chosen the next fiber to break, we remove
            # this tension and calculate it again
            σ[i] = 1
        
        # else
        # The fiber is BROKEN and should stay BROKEN, 
        # so there is no need to do anything
        end
    end
end

function findNextFiber(σ, x)
    # σ is the relative tension of the fiber if x had been 1
    # If σ of a fiber is 2, this just means that it is under
    # twice as much tension as a fiber of σ=1. But in order to
    # find what fiber will break first, we need to take into
    # account how much tension the fiber can handle. We do this
    # by dividing by x.
    return @fastmath argmax(σ ./ x)
end

function break_fiber(i::Int64, status::Vector{Int64}, σ::Vector{Float64})
    @inbounds begin
        status[i] = 0#BROKEN
        σ[i] = 0
    end
end

function update_σ(status::Vector{Int64}, σ::Vector{Float64},
    neighbours::Array{Int64, 2},
    neighbourhoods::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_dimensions::Vector{Int64},
    rel_pos_x::Vector{Int64},
    rel_pos_y::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64};
    neighbourhood_rule::String="UNR")
    # Explores the plane, identifies all the clusters, their sizes
    # and outlines

    # Cluster id
    c::Int64 = 0
    # The id of a spanning cluster. If there are none, it is -1
    spanning_cluster::Int64 = -1

    # For every fiber in the plane
    @fastmath @inbounds for i in eachindex(status)
        # If it is broken and unexplored
        if status[i] == 0#BROKEN
            # We have found a new cluster!
            # increase number of clusters
            c += 1
            # assign fiber i to this new cluster
            status[i] = c
            # explore the new cluster
            spanning_cluster = explore_cluster_at(i, c, status, neighbours, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored)
            # We should now have updated cluster_outline,
            # and with that we can update sigma for one cluster
            update_cluster_outline_stress(c,status,σ, cluster_size, cluster_outline, cluster_outline_length, neighbourhoods, neighbourhood_rule)
        end
    end
    spanning_cluster_size = spanning_cluster==-1 ? -1 : cluster_size[spanning_cluster]
    return c, spanning_cluster, spanning_cluster_size
end

function explore_cluster_at(i::Int64, c::Int64,
    status::Vector{Int64},
    neighbours::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_dimensions::Vector{Int64},
    rel_pos_x::Vector{Int64},
    rel_pos_y::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64})
    # We explore the cluster of broken fibers and
    # map the border of the cluster
    @fastmath begin
        # Number of fibers in current cluster that has been explored
        nr_explored::Int64 = 0
        # Number of unexplored fibers that we know of.
        nr_unexplored::Int64 = 1 # Starts at one because of i
        # These are the actual unexplored. We just overwrite new values
        unexplored[1] = i
        # This value will be continuesly updated as we explore the cluster
        cluster_size[c] = 1
        # This value will be continuesly updated as we explore the cluster
        cluster_outline_length[c] = 0

        # Cluster dimensions
        # [max_x, min_x, max_y, min_y]
        fill!(cluster_dimensions, 0) # Reset cluster dimensions

        # While there are still unexplored fibers in the cluster
        while nr_unexplored > nr_explored
            # Preemptively count this fiber as explored (because 1 indexing)
            nr_explored += 1
            # Get the next unexplored fiber
            current_fiber = unexplored[nr_explored]
            # Go through all neighbours of the fiber
            nr_unexplored = check_neighbours(current_fiber, nr_unexplored, c, status, neighbours, cluster_size,
                                             cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline,
                                             cluster_outline_length, unexplored)
        end

        L::Int64 = isqrt(length(status)) 
        if spanning(L, cluster_dimensions)
            return c
        else 
            return -1
        end

    end
end

function spanning(L::Int64, cluster_dimensions::Vector{Int64})
    # Cluster dimensions
    # [max_x, min_x, max_y, min_y]
    # L-1 because the relative coordinates in the cluster start at 0,0
    # NB! Once the cluster is spanning, the dimension is no longer reliable because of
    # periodicity.
    return cluster_dimensions[1] - cluster_dimensions[2] >= L-1 || cluster_dimensions[3] - cluster_dimensions[4] >= L-1
end

function check_neighbours(current_fiber::Int64, nr_unexplored::Int64, c::Int64,
    status::Vector{Int64},
    neighbours::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_dimensions::Vector{Int64},
    rel_pos_x::Vector{Int64},
    rel_pos_y::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64})
    # Here we check the neighbour fibers of a given fiber
    # We want to see if the neighbours are broken, making them
    # part of the cluster, or if they are alive, in which case
    # we need to add them to the border of the cluster.

    @fastmath for (i, neighbour_fiber) in enumerate(neighbours[current_fiber, :])
        # Status of neighbour fiber
        s::Int64 = status[neighbour_fiber]
        # If this adjecent fiber is is BROKEN,
        if s == 0#BROKEN
            # We then have to add this to the list of unexplored 
            nr_unexplored = add_unexplored(neighbour_fiber, unexplored, nr_unexplored)
            # We set this fiber to be part of the current cluster
            status[neighbour_fiber] = c
            # increase the cluster size
            cluster_size[c] += 1
            # and increase the cluster dimensions depending on what direction we explored In
            store_possition(current_fiber, neighbour_fiber, i, cluster_dimensions, rel_pos_x, rel_pos_y)
        # In some situations, a fiber will be part of the border of
        # two different clusters, so we check for ALIVE or PAST_BORDER
        elseif s == -1 || s == -3 #ALIVE || PAST_BORDER
            # We have to change to CURRENT_BORDER so that
            # we don't count it again since a fiber will often
            # be checked multiple times
            status[neighbour_fiber] = -2 #CURRENT_BORDER
            cluster_outline_length[c] += 1
            cluster_outline[cluster_outline_length[c]] = neighbour_fiber
        end
    end
    return nr_unexplored
end

function store_possition(current_fiber::Int64, neighbour_fiber::Int64, 
    direction::Int64, 
    cluster_dimensions::Vector{Int64},
    rel_pos_x::Vector{Int64},
    rel_pos_y::Vector{Int64})
    # cluster_dimensions = [max_x, min_x, max_y, min_y]
    movement = [1,-1,1,-1] # direction

    pos = direction<3 ? rel_pos_x : rel_pos_y # If direction is 1 or 2, it is the x direction we use

    pos[neighbour_fiber] = pos[current_fiber] + movement[direction]
    if abs(pos[neighbour_fiber]) > abs(cluster_dimensions[direction])
        cluster_dimensions[direction] = pos[neighbour_fiber]
    end
end

function add_unexplored(i::Int64, unexplored::Vector{Int64}, nr_unexplored::Int64)
    @fastmath @inbounds begin 
        nr_unexplored += 1
        unexplored[nr_unexplored] = i
        return nr_unexplored
    end
end

function apply_to_neighbourhood(f::Function,
    status::Vector{Int64},
    cluster_outline::Vector{Int64},
    neighbourhoods::Array{Int64, 2})
    # For every fiber in the cluster outline, take the 3x3 matrix around the fiber and 
    # put it into the function f
    value::Vector{Int64} = zeros(Int64, length(cluster_outline))
    @fastmath @inbounds @simd for i in 1:length(cluster_outline)
        outline_fiber = cluster_outline[i]
        value[i] = f(status[neighbourhoods[outline_fiber, :]])
    end
    return value
end

function alive_fibers_in_neighbourhood(m::Vector{Int64})
    # Take a 3x3 matrix around a fiber and count how many are alive
    alive_fibers = 1 # We count the center fiber as well because the math doesn't like 0
    @fastmath @inbounds @simd for f in m
        if f<0
            alive_fibers +=1
        end
    end
    return alive_fibers
end

function update_cluster_outline_stress(c::Int64,
    status::Vector{Int64},
    σ::Vector{Float64},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    neighbourhoods::Array{Int64, 2},
    neighbourhood_rule::String)
    # Apply the appropreate amount of stress to the fibers
    if neighbourhood_rule == "UNR"
        # With the Uniform neighbourhood rule, we can apply a simple stress
        apply_simple_stress(c, status, σ, cluster_size, cluster_outline, cluster_outline_length)
        return
    else
        # But with more complex rules, we need to do it in two steps
        # First a calculation to find the fiber strengths (As a function of their neighbourhood), and then apply the stress
        fiber_strengths = ones(Float64)
        if neighbourhood_rule == "CNR"
            outline_neihbourhoods = apply_to_neighbourhood(neighbourhoodToInt, status, cluster_outline[1:cluster_outline_length[c]], neighbourhoods)
            fiber_strengths = neighbourhoodStrengths[outline_neihbourhoods]
        elseif neighbourhood_rule == "SNR"
            fiber_strengths = apply_to_neighbourhood(alive_fibers_in_neighbourhood, status, cluster_outline[1:cluster_outline_length[c]], neighbourhoods)
        else
            error("Unknown neighbourhood rule")
        end

        apply_stress(c, status, σ, cluster_size, cluster_outline, cluster_outline_length, fiber_strengths)
    end
end

function apply_simple_stress(c::Int64,
    status::Vector{Int64},
    σ::Vector{Float64},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64}
    )
    added_stress = cluster_size[c]/cluster_outline_length[c]
    @inbounds @simd for i in 1:cluster_outline_length[c]
        fiber = cluster_outline[i]
        σ[fiber] += added_stress
        status[fiber] = -3 #PAST_BORDER
    end
end

function apply_stress(c::Int64,
    status::Vector{Int64},
    σ::Vector{Float64},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    fiber_strengths::Vector{Float64}
    )
    # See page 26 in Jonas Tøgersen Kjellstadli's doctoral theses, 2019:368
    α = 2.0 # High alpha means that having neighbours is more important
    C = 1 / sum(fiber_strengths .^(-α+1)) # Normalization constant
    g = C .* fiber_strengths .^(-α) # Normalization factor
    @simd for i in 1:cluster_outline_length[c]
        fiber = cluster_outline[i]
        added_stress =  cluster_size[c]*fiber_strengths[i]*g[i]
        σ[fiber] += added_stress
        status[fiber] = -3 #PAST_BORDER
    end
end