using BenchmarkTools
using Random
using ProgressMeter

"""
Simulate breakage of an LxL fiber matrix
"""

# Define constants
ALIVE::Int64 = -1 # A fiber that has not yet been broken
CURRENT_BORDER::Int64 = -2 # The border of the current cluster being explored
PAST_BORDER::Int64 = -3 # The border of a cluster that has been explored
BROKEN::Int64 = 0 # A broken fiber that has not been explored
# Status larger than 0 represents the cluster id that the fiber belongs to

function fillAdjacent(L::Int64)
    N = L*L
    adjacent = zeros(Int64, N, 4)
    for i in 1:N
        adjacent[i, :] = findAdjacent(i, L)
    end
    return adjacent
end

function findAdjacent(i::Int64, L::Int64)
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
    a = [
        i-1, # Up
        i+1, # Down
        i+L, # Right
        i-L  # Left
        ]

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

    return a 
end

function resetClusters(status::Vector{Int64}, σ::Vector{Float64})
    # After having explored and assigned numbers to 
    # all the fibers indicating their state, we now 
    # want to reset them to either be BROKEN or ALIVE
    for i in eachindex(status)
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
    return argmax(σ ./ x)
end

function break_fiber(i::Int64, status::Vector{Int64}, σ::Vector{Float64})
    status[i] = BROKEN
    σ[i] = 0
end

function update_σ(status::Vector{Int64}, σ::Vector{Float64},
    adjacent::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64})
    # Explores the plane, identifies all the clusters, their sizes
    # and outlines

    # Cluster id
    c::Int64 = 0


    # For every fiber in the plane
    for i in eachindex(status)
        # If it is broken and unexplored
        if status[i] == BROKEN
            # We have found a new cluster!
            # increase number of clusters
            c += 1
            # assign fiber i to this new cluster
            status[i] = c
            # explore the new cluster
            explore_cluster_at(i, c, status, adjacent, cluster_size, cluster_outline, cluster_outline_length, unexplored)
            # We should now have updated cluster_outline,
            # and with that we can update sigma for one cluster
            update_cluster_outline_stress(c,status,σ, cluster_size, cluster_outline, cluster_outline_length)
        end
    end
end

function explore_cluster_at(i::Int64, c::Int64,
    status::Vector{Int64},
    adjacent::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64})
    # We explore the cluster of broken fibers and
    # map the border of the cluster

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
        nr_unexplored = check_neighbours(current_fiber,nr_unexplored, c, status, adjacent, cluster_size, cluster_outline, cluster_outline_length, unexplored)
    end
end

function check_neighbours(current_fiber::Int64, nr_unexplored::Int64, c::Int64,
    status::Vector{Int64},
    adjacent::Array{Int64, 2},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64},
    unexplored::Vector{Int64})
    # Here we check the adjacent fibers of a given fiber
    # We want to see if the neighbours are broken, making them
    # part of the cluster, or if they are alive, in which case
    # we need to add them to the border of the cluster.

    for adj_fiber in adjacent[current_fiber, :]
        # Status of adjacent fiber
        s::Int64 = status[adj_fiber]
        # If this adjecent fiber is is BROKEN,
        if s == BROKEN
            # We then have to add this to the list of unexplored 
            nr_unexplored = add_unexplored(adj_fiber, unexplored, nr_unexplored)
            # We set this fiber to be part of the current cluster
            status[adj_fiber] = c
            # and increase the cluster size
            cluster_size[c] += 1
        # In some situations, a fiber will be part of the border of
        # two different clusters, so we check for ALIVE or PAST_BORDER
        elseif s == ALIVE || s == PAST_BORDER
            # We have to change to CURRENT_BORDER so that
            # we don't count it again since a fiber will often
            # be checked multiple times
            status[adj_fiber] = CURRENT_BORDER
            cluster_outline_length[c] += 1
            cluster_outline[cluster_outline_length[c]] = adj_fiber
        end
    end
    return nr_unexplored
end

function add_unexplored(i::Int64, unexplored::Vector{Int64}, nr_unexplored::Int64)
    nr_unexplored += 1
    unexplored[nr_unexplored] = i
    return nr_unexplored
end

function update_cluster_outline_stress(c::Int64,
    status::Vector{Int64},
    σ::Vector{Float64},
    cluster_size::Vector{Int64},
    cluster_outline::Vector{Int64},
    cluster_outline_length::Vector{Int64})

    for i in 1:cluster_outline_length[c]
        fiber = cluster_outline[i]
        added_stress = cluster_size[c]/cluster_outline_length[c]

        σ[fiber] += added_stress
        status[fiber] = PAST_BORDER
    end
end

function main(L)
    N = L*L # Number of fibers
    x = rand(Float64, N) # Max extension
    σ  = ones(Float64, N) # Relative tension
    adjacent = fillAdjacent(L)
    status = fill(-1, N)
    cluster_size = zeros(Int64, N)
    cluster_outline = zeros(Int64, N)
    cluster_outline_length = zeros(Int64, N)
    unexplored = zeros(Int64, N)

    # Do what you want
    @showprogress for step in 1:N
        i = findNextFiber(σ, x)
        resetClusters(status, σ)
        break_fiber(i, status, σ)
        update_σ(status,σ,adjacent, cluster_size, cluster_outline, cluster_outline_length, unexplored)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main(128)
end