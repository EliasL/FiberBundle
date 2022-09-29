using BenchmarkTools
using Random
using ProgressMeter
using SparseArrays

include("NeighbourhoodWeighting.jl")

"""
Simulate breakage of an LxL fiber matrix
"""

# Define constants
ALIVE::Int64 = -1 # A fiber that has not yet been broken
CURRENT_BORDER::Int64 = -2 # The border of the current cluster being explored
PAST_BORDER::Int64 = -3 # The border of a cluster that has been explored
BROKEN::Int64 = 0 # A broken fiber that has not been explored
# Status larger than 0 represents the cluster id that the fiber belongs to

 #
#x#
 #
NEIGHBOURS(i, L) = [
    i-1, # Up
    i+1, # Down
    i+L, # Right
    i-L  # Left
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
    i,     # Itself
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
    # The structure of the matrix will be
    # 1 4 7
    # 2 5 8
    # 3 6 9
    # If the checks are confusing, use the gid above to
    # see what happens.
    # mod1(7,3) = mod1(1,3) = 1
    # mod1(3,3) = mod1(9,3) = 3
    N = L*L
    # Fill in adjacent
    a = adjacent

    if length(a) == 4
        # Up check
        if mod1(i, L) == 1
            a[1] += L
        end
        # Down check
        if mod1(i, L) == L
            a[2] -= L
        end
        # Right check
        if i>N-L
            a[3] -= N
        end
        # Left check
        if i<=L
            a[4] += N
        end
    elseif length(a) == 9
        # Up check
        if mod1(i, L) == 1
            a[1] += L
            a[4] += L
            a[7] += L
        end
        # Down check
        if mod1(i, L) == L
            a[3] -= L
            a[6] -= L
            a[9] -= L
        end
        # Right check
        if i>N-L
            a[7] -= N
            a[8] -= N
            a[9] -= N
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
        if s > BROKEN
            # In this case, we know that it is broken
            # and should be reset
            status[i] = BROKEN

        elseif s < BROKEN
            # This means that this fiber was set as either
            # ALIVE, CURRENT_BORDER or PAST_BORDER.
            # Now we reset it to being just ALIVE
            status[i] = ALIVE
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
        status[i] = BROKEN
        σ[i] = 0
    end
end

function update_σ(status::Vector{Int64}, σ::Vector{Float64},
    neighbours::Array{Int64, 2},
    neighbourhoods::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64};
    use_neighbourhood_rules::Bool=false)
    # Explores the plane, identifies all the clusters, their sizes
    # and outlines

    # Cluster id
    c::Int64 = 0


    # For every fiber in the plane
    @fastmath @inbounds for i in eachindex(status)
        # If it is broken and unexplored
        if status[i] == BROKEN
            # We have found a new cluster!
            # increase number of clusters
            c += 1
            # assign fiber i to this new cluster
            status[i] = c
            # explore the new cluster
            explore_cluster_at(i, c, status, neighbours, cluster_size, cluster_outline, cluster_outline_length, unexplored)
            # We should now have updated cluster_outline,
            # and with that we can update sigma for one cluster
            if use_neighbourhood_rules
                update_cluster_otline_stress_with_neighbourhood_rules(c,status,σ, cluster_size, cluster_outline, cluster_outline_length, neighbourhoods)
            else
                update_cluster_outline_stress(c,status,σ, cluster_size, cluster_outline, cluster_outline_length)
            end
        end
    end
    return c
end

function explore_cluster_at(i::Int64, c::Int64,
    status::Vector{Int64},
    neighbours::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64})
    # We explore the cluster of broken fibers and
    # map the border of the cluster
    @fastmath @inbounds begin
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

        # While there are still unexplored fibers in the cluster
        while nr_unexplored > nr_explored
            # Preemptively count this fiber as explored (because 1 indexing)
            nr_explored += 1
            # Get the next unexplored fiber
            current_fiber = unexplored[nr_explored]
            # Go through all neighbours of the fiber
            nr_unexplored = check_neighbours(current_fiber,nr_unexplored, c, status, neighbours, cluster_size, cluster_outline, cluster_outline_length, unexplored)
        end
    end
end

function check_neighbours(current_fiber::Int64, nr_unexplored::Int64, c::Int64,
    status::Vector{Int64},
    neighbours::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64})
    # Here we check the neighbour fibers of a given fiber
    # We want to see if the neighbours are broken, making them
    # part of the cluster, or if they are alive, in which case
    # we need to add them to the border of the cluster.

    @fastmath @inbounds @simd for neighbour_fiber in neighbours[current_fiber, :]
        # Status of neighbour fiber
        s::Int64 = status[neighbour_fiber]
        # If this adjecent fiber is is BROKEN,
        if s == BROKEN
            # We then have to add this to the list of unexplored 
            nr_unexplored = add_unexplored(neighbour_fiber, unexplored, nr_unexplored)
            # We set this fiber to be part of the current cluster
            status[neighbour_fiber] = c
            # and increase the cluster size
            cluster_size[c] += 1
        # In some situations, a fiber will be part of the border of
        # two different clusters, so we check for ALIVE or PAST_BORDER
        elseif s == ALIVE || s == PAST_BORDER
            # We have to change to CURRENT_BORDER so that
            # we don't count it again since a fiber will often
            # be checked multiple times
            status[neighbour_fiber] = CURRENT_BORDER
            cluster_outline_length[c] += 1
            cluster_outline[cluster_outline_length[c]] = neighbour_fiber
        end
    end
    return nr_unexplored
end

function add_unexplored(i::Int64, unexplored::Vector{Int64}, nr_unexplored::Int64)
    @fastmath @inbounds begin 
        nr_unexplored += 1
        unexplored[nr_unexplored] = i
        return nr_unexplored
    end
end

function update_cluster_outline_stress(c::Int64,
    status::Vector{Int64},
    σ::Vector{Float64},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64})

    @fastmath @inbounds @simd for i in 1:cluster_outline_length[c]
        fiber = cluster_outline[i]
        added_stress = cluster_size[c]/cluster_outline_length[c]

        σ[fiber] += added_stress
        status[fiber] = PAST_BORDER
    end
end

function get_id_of_neighbourhoods_of_outline(
    status::Vector{Int64},
    cluster_outline::Vector{Int64},
    neighbourhoods::Array{Int64, 2})
    return map(o -> neighbourhoodToInt(status[neighbourhoods[o, :]]), cluster_outline)
end

function update_cluster_otline_stress_with_neighbourhood_rules(c::Int64,
    status::Vector{Int64},
    σ::Vector{Float64},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    neighbourhoods::Array{Int64, 2})

    # Identifies the id of the neighbourhood of each outline fiber
    outline_neihbourhoods = get_id_of_neighbourhoods_of_outline(status, cluster_outline[1:cluster_outline_length[c]], neighbourhoods)
    outline_strengths = neighbourhoodStrengths[outline_neihbourhoods]

    # initial_stress is an interesting parameter. If max strength is 11,
    # then fibers whos neighbourhood has a missing corner will experience
    # no stress at all. If the initial_stress is larger, for example 100,
    # it means that the difference between the neighbourhood configurations
    # is less significant. As initial_stress -> inf, the model will approach
    # the LLS
    initial_stress = 12
    outline_stress = initial_stress .- outline_strengths
    total_outline_stress = sum(outline_stress)
    average_outline_stress = total_outline_stress/cluster_outline_length[c]

    for i in 1:cluster_outline_length[c]
        fiber = cluster_outline[i]
        added_stress = cluster_size[c] / cluster_outline_length[c] * outline_stress[i] / average_outline_stress
        σ[fiber] += added_stress
        status[fiber] = PAST_BORDER
    end
end