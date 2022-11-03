using BenchmarkTools
using Random
using ProgressMeter
using StaticArrays

include("support/neighbourhoodWeighting.jl")
include("support/distributions.jl")

"""
Simulate breakage of an LxL fiber matrix
"""

# Fiber bundle

#An idea is to store variables in a mutable vector so that the struct can stay static (Not mulable)
Base.@kwdef mutable struct FB{l, n, F<:AbstractFloat, I<:Integer}
    L::I = l
    N::I = n
    α::F = 2
    nr::String = "UNR"
    x::MVector{n,F} = MVector{n}(zeros(F, n))
    neighbours::SMatrix{n, 4, I} = SMatrix{n, 4}(fillAdjacent(l, NEIGHBOURS))
    neighbourhoods::SMatrix{n, 8, I} = SMatrix{n, 8}(fillAdjacent(l, NEIGHBOURHOOD))
    movement::SVector{4,I} = SVector{4}([1,-1,1,-1]) # this is dependent on the order of neighbours...
    current_neighbourhood::MVector{8, I} = MVector{8}(zeros(I, 8))
    neighbourhood_values::MVector{n,I} = MVector{n}(zeros(I, n))

    # These values are reset for each step
    σ::MVector{n,F} = MVector{n}(ones(F, n)) # Relative tension
    tension::MVector{n,F} = MVector{n}(zeros(F, n))
    max_σ::F = 0.0
    status::MVector{n,I} = MVector{n}(fill(I(-1), n))
    c::I=0
    spanning_cluster_id::I=-1
    cluster_size::MVector{n,I} = MVector{n}(zeros(I, n))
    cluster_outline_length::MVector{n,I} = MVector{n}(zeros(I, n))
    # These values are reset for each cluster
    cluster_outline::MVector{n,I} = MVector{n}(zeros(I, n))
    unexplored::MVector{n,I} = MVector{n}(zeros(I, n))
    # Relative possition of every fiber with respect to it's cluster
    rel_pos_x::MVector{n,I} = MVector{n}(zeros(I, n))
    rel_pos_y::MVector{n,I} = MVector{n}(zeros(I, n))
    cluster_dimensions::MVector{4,I} = MVector{4}(zeros(I, 4))
end

# Fiber bundle storage
Base.@kwdef struct FBS{division, n, F<:AbstractFloat, I<:Integer}
    # These arrays store one value for each step
    most_stressed_fiber::MVector{n,F} = MVector{n}(zeros(F, n))
    nr_clusters::MVector{n,I} = MVector{n}(zeros(I, n))
    largest_cluster::MVector{n,I} = MVector{n}(zeros(I, n))
    largest_perimiter::MVector{n,I} = MVector{n}(zeros(I, n))

    # We want to store some samples of the processed
    # I'm thinking at 10%, 20%, ... 90% done would work
    # ie, 9 images
    status_storage = zeros(I, division-1, n)
    tension_storage = zeros(F, division-1, n)

    spanning_cluster_state_storage = zeros(I, n)
    spanning_cluster_tension_storage = zeros(I, n)
    spanning_cluster_size_storage = 0
    spanning_cluster_perimiter_storage = 0
    spanning_cluster_has_not_been_found = true
    spanning_cluster_step = 0
    # If N=100 Steps to store is now [90, 80, ... , 10]
    steps_to_store = [round(I,n/division * i) for i in 1:division-1]
    storage_index = 1
end

function get_fb(L; α=2, t=0, nr="UNR", dist="Uniform", without_storage=false)
    N=L*L
    x = nothing
    if dist == "Uniform"
        distribution_function = get_uniform_distribution(t)
        x = distribution_function(N)
    elseif isa(dist, Function)
        x = dist(N)
    elseif isa(dist, AbstractVector)
        x = dist
    else
        error("No distribution found! Got: $dist")
    end

    fb = FB{L, N, Float64, Int64}(α=α, nr=nr, x=x)
    if without_storage
        return fb
    else
        return fb, FBS{10, N, Float64, Int64}()
    end
end

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

function resetBundle!(b::FB)
    #NB This does not completely reset the cluster, only
    # values that should be reset every step
    
    #Reset cluster id/number 
    b.c = 0
    # The id of a spanning cluster. If there are none, it is -1
    b.spanning_cluster_id = -1

    resetClusters!(b)

    reset_relative_possition!(b)
end

function resetClusters!(b::FB)
    # After having explored and assigned numbers to 
    # all the fibers indicating their state, we now 
    # want to reset them to either be BROKEN or ALIVE

    for i in eachindex(b.status)
        s = b.status[i]

        # The statuses of the fibers will indicate
        # what cluster they belong to. If this fiber
        # has had a status larger than 0, it means that
        # it has belonged to a cluster. 
        if s > 0#BROKEN
            # In this case, we know that it is broken
            # and should be reset
            b.status[i] = 0#BROKEN

        elseif s < 0#BROKEN
            # This means that this fiber was set as either
            # ALIVE, CURRENT_BORDER or PAST_BORDER.
            # Now we reset it to being just ALIVE
            b.status[i] = -1 #ALIVE
            # We have already added stress to the boarder, so now
            # AFTER having chosen the next fiber to break, we remove
            # this tension and calculate it again
            b.σ[i] = 1
        
        # else
        # The fiber is BROKEN and should stay BROKEN, 
        # so there is no need to do anything
        end
    end
end

function update_tension!(b::FB)# σ is the relative tension of the fiber if x had been 1
    # If σ of a fiber is 2, this just means that it is under
    # twice as much tension as a fiber of σ=1. But in order to
    # find what fiber will break first, we need to take into
    # account how much tension the fiber can handle. We do this
    # by dividing by x.
    for i in eachindex(b.σ)
        b.tension[i] = b.σ[i] / b.x[i]
    end
end

function findNextFiber!(b::FB)
    update_tension!(b)
    i = argmax(b.tension)
    b.max_σ = b.tension[i]
    return i
end

function break_fiber!(i::Int, b::FB)
    b.status[i] = 0#BROKEN
    b.σ[i] = 0
end

function reset_relative_possition!(b::FB)
    # Note about the relative possition:
    # If a fiber boarders two clusters, it will have the relative
    # possition of the cluster explored last.
    fill!(b.rel_pos_x, 0)
    fill!(b.rel_pos_y, 0)
end

function update_σ!(b::FB)
    # Explores the plane, identifies all the clusters, their sizes
    # and outlines


    # For every fiber in the plane
    for i in eachindex(b.status)
        #@logmsg σUpdateLog "Breaking fiber $i"
        # If it is broken and unexplored
        if b.status[i] == 0#BROKEN
            # We have found a new cluster!
            # increase number of clusters
            b.c += 1
            # assign fiber i to this new cluster
            b.status[i] = b.c
            # explore the new cluster
            #@logmsg σUpdateLog "Exploring cluster"
            found_spanning_cluster = explore_cluster_at!(i, b)
            if found_spanning_cluster && b.spanning_cluster_id == -1
                b.spanning_cluster_id = b.c
            end
            # We should now have updated cluster_outline,
            # and with that we can update sigma for one cluster
            
            update_cluster_outline_stress!(b)
        end
    end
end

function reset_cluster_dimensions!(b::FB)

    # Cluster dimensions
    # [max_x, min_x, max_y, min_y]
    fill!(b.cluster_dimensions, 0) # Reset cluster dimensions
end

function explore_cluster_at!(i::Int64, b::FB)
    # We explore the cluster of broken fibers and
    # map the border of the cluster
    # Number of fibers in current cluster that has been explored
    nr_explored::Int64 = 0
    # Number of unexplored fibers that we know of.
    nr_unexplored::Int64 = 1 # Starts at one because of i
    # These are the actual unexplored. We just overwrite new values
    b.unexplored[1] = i
    # This value will be continuesly updated as we explore the cluster
    b.cluster_size[b.c] = 1
    # This value will be continuesly updated as we explore the cluster
    b.cluster_outline_length[b.c] = 0

    reset_cluster_dimensions!(b)

    # While there are still unexplored fibers in the cluster
    while nr_unexplored > nr_explored

        #@logmsg clusterLog "Explored $nr_unexplored / $nr_explored"
        # Preemptively count this fiber as explored (because 1 indexing)
        nr_explored += 1
        # Get the next unexplored fiber
        current_fiber = b.unexplored[nr_explored]
        # Go through all neighbours of the fiber
        nr_unexplored = check_neighbours!(current_fiber, nr_unexplored, b)
    end
    
    if spanning(b)
        return true
    else 
        return false
    end


end

function spanning(b::FB)
    # Cluster dimensions
    # [max_x, min_x, max_y, min_y]
    # L-1 because the relative coordinates in the cluster start at 0,0
    # NB! Once the cluster is spanning, the dimension is no longer reliable because of
    # periodicity.
    return b.cluster_dimensions[1] - b.cluster_dimensions[2] >= b.L-1 ||
           b.cluster_dimensions[3] - b.cluster_dimensions[4] >= b.L-1
end

function check_neighbours!(current_fiber::Int64, nr_unexplored::Int64, b::FB)
    # Here we check the neighbour fibers of a given fiber
    # We want to see if the neighbours are broken, making them
    # part of the cluster, or if they are alive, in which case
    # we need to add them to the border of the cluster.
    for (i, neighbour_fiber) in enumerate(view(b.neighbours, current_fiber, :))
        # Status of neighbour fiber
        s::Int64 = b.status[neighbour_fiber]
        # If this adjecent fiber is is BROKEN,
        if s == 0#BROKEN
            # We then have to add this to the list of unexplored 
            nr_unexplored = add_unexplored!(neighbour_fiber, nr_unexplored, b)
            # We set this fiber to be part of the current cluster
            b.status[neighbour_fiber] = b.c
            # increase the cluster size
            b.cluster_size[b.c] += 1
            # and increase the cluster dimensions depending on what direction we explored In
            store_possition!(current_fiber, neighbour_fiber, i, b)
        # In some situations, a fiber will be part of the border of
        # two different clusters, so we check for ALIVE or PAST_BORDER
        elseif s == -1 || s == -3 #ALIVE || PAST_BORDER
            # We have to change to CURRENT_BORDER so that
            # we don't count it again since a fiber will often
            # be checked multiple times
            b.status[neighbour_fiber] = -2 #CURRENT_BORDER
            b.cluster_outline_length[b.c] += 1
            b.cluster_outline[b.cluster_outline_length[b.c]] = neighbour_fiber
        end
    end
    return nr_unexplored
end


function store_possition!(current_fiber::Int64, neighbour_fiber::Int64, 
    direction::Int64, 
    b::FB)
    # cluster_dimensions = [max_x, min_x, max_y, min_y]

    # Copy over the possition to the neighbour
    b.rel_pos_x[neighbour_fiber] = b.rel_pos_x[current_fiber]
    b.rel_pos_y[neighbour_fiber] = b.rel_pos_y[current_fiber]

    #Then add the movement to the neighbour
    xOrY = direction<=2 ? b.rel_pos_x : b.rel_pos_y # If direction is 1 or 2, it is the x direction we use
    xOrY[neighbour_fiber] += b.movement[direction]

    # If this is a new max, then we save it
    if xOrY[neighbour_fiber]*b.movement[direction] > b.cluster_dimensions[direction]*b.movement[direction]
        b.cluster_dimensions[direction] = xOrY[neighbour_fiber]
    end
end

function add_unexplored!(i::Int64, nr_unexplored::Int64, b::FB)
    nr_unexplored += 1
    b.unexplored[nr_unexplored] = i
    return nr_unexplored
end

function apply_to_neighbourhood!(f::Function, b::FB)
    # For every fiber in the cluster outline, take the 3x3 matrix around the fiber and 
    # put it into the function f

    # At the time this function is run, b.c is the current cluster
    for i in 1:b.cluster_outline_length[b.c]
        update_current_neighbourhood!(i, b)
        b.neighbourhood_values[i] = f(b.current_neighbourhood)
    end
end

function update_current_neighbourhood!(i::Int64, b::FB)
    for j in 1:8
        b.current_neighbourhood[j] = b.status[b.neighbourhoods[b.cluster_outline[i], j]]
    end
end


function alive_fibers_in_neighbourhood(n::AbstractVector{Int64})
    # Take a neighbourhood and count the number of alive fibers
    alive_fibers = 1 # We count the center fiber as well because the math doesn't like 0
    for f in n
        if f<0
            alive_fibers +=1
        end
    end
    return alive_fibers
end

function update_cluster_outline_stress!(b::FB)
    # Apply the appropreate amount of stress to the fibers
    if b.nr == "UNR"
        # With the Uniform neighbourhood rule, we can apply a simple stress
        apply_simple_stress(b)
        return
    else
        # But with more complex rules, we need to do it in two steps
        # First a calculation to find the fiber strengths (As a function of their neighbourhood), and then apply the stress
        if b.nr == "CNR"
            apply_to_neighbourhood!(neighbourhoodToStrength, b)
        elseif b.nr == "SNR"
            apply_to_neighbourhood!(alive_fibers_in_neighbourhood, b)
        else
            #@debug "Unknown neighbourhood rule: $neighbourhood_rule"
            error("Unknown neighbourhood rule")
        end

        apply_stress(b)
    end
end

function apply_simple_stress(b::FB)
    added_stress = b.cluster_size[b.c]/b.cluster_outline_length[b.c]
    for i in 1:b.cluster_outline_length[b.c]
        fiber = b.cluster_outline[i]
        b.σ[fiber] += added_stress
        b.status[fiber] = -3 #PAST_BORDER
    end
end


function apply_stress(b::FB)
    # See page 26 in Jonas Tøgersen Kjellstadli's doctoral theses, 2019:368
    # High alpha means that having neighbours is more important
    C=0.0
    for i in 1:b.cluster_outline_length[b.c]
        C += b.neighbourhood_values[i] ^(-b.α+1)
    end
    C = 1/C # A normalization constant
    for i in 1:b.cluster_outline_length[b.c]
        fiber = b.cluster_outline[i]
        g = C * b.neighbourhood_values[i] ^(-b.α)
        added_stress =  b.cluster_size[b.c]*b.neighbourhood_values[i]*g
        b.σ[fiber] += added_stress
        b.status[fiber] = -3 #PAST_BORDER
    end
end

#fb, storage = get_fb(8)