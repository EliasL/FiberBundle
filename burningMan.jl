using BenchmarkTools
using Random
using ProgressMeter
using StaticArrays
using Profile
#using PProf

include("support/neighbourhoodWeighting.jl")
include("support/distributions.jl")

"""
Simulate breakage of an LxL fiber matrix
"""

# Fiber bundle

#An idea is to store variables in a mutable vector so that the struct can stay static (Not mulable)
Base.@kwdef mutable struct FB{F<:AbstractFloat, I<:Integer}
    L::I
    N::I = L*L
    α::F = 2
    nr::String = "LLS"
    x::Vector{F} = zeros(F, N)
    neighbours::Matrix{I} = fillAdjacent(L, NEIGHBOURS)
    neighbourhoods::Matrix{I} = fillAdjacent(L, NEIGHBOURHOOD)
    movement::SVector{4,I} = SVector{4}([1,-1,-1,1]) # this is dependent on the order of neighbours...
    current_neighbourhood::Vector{I} = zeros(I, 8)
    neighbourhood_values::Vector{I} = zeros(I, N)

    # These values are reset for each step
    σ::Vector{F} = ones(F, N) # Relative tension
    tension::Vector{F} = zeros(F, N)
    max_σ::F = 0.0
    status::Vector{I} = fill(I(-1), N)
    current_step::I = 0
    break_sequence::Vector{I} = zeros(I, N)
    c::I = 0 # nr_clusters / current_cluster_id
    spanning_cluster_id::I = -1
    cluster_size::Vector{I} = zeros(I, N)
    cluster_outline_length::Vector{I} = zeros(I, N)
    # These values are reset for each cluster
    cluster_outline::Vector{I} = zeros(I, N)
    cluster_cm_x::Vector{F} = zeros(F, N)
    cluster_cm_y::Vector{F} = zeros(F, N)
    unexplored::Vector{I} = zeros(I, N)
    # Relative possition of every fiber with respect to it's cluster
    rel_pos_x::Vector{I} = zeros(I, N)#Vector{I}(undef, N)
    rel_pos_y::Vector{I} = zeros(I, N)#Vector{I}(undef, N)
    cluster_dimensions::Vector{I} = zeros(I, 4)
end

# Fiber bundle storage
Base.@kwdef mutable struct FBS{F<:AbstractFloat, I<:Integer}
    division::I
    N::I
    # These arrays store one value for each step
    most_stressed_fiber::Vector{F} = zeros(F, N)
    nr_clusters::Vector{I} = zeros(I, N)
    largest_cluster::Vector{I} = zeros(I, N)
    largest_perimiter::Vector{I} = zeros(I, N)

    # We want to store some samples of the processed
    # I'm thinking at 10%, 20%, ... 90% done would work
    # ie, 9 images
    status_storage = zeros(I, division-1, N)
    tension_storage = zeros(F, division-1, N)

    spanning_cluster_state_storage = zeros(I, N)
    spanning_cluster_tension_storage = zeros(F, N)
    spanning_cluster_size_storage::I = 0
    spanning_cluster_perimiter_storage::I = 0
    spanning_cluster_has_been_found::Bool = false
    spanning_cluster_step::I = 0
    # If N=100 Steps to store is now [90, 80, ... , 10]
    steps_to_store::Vector{I} = [round(I,N/division * i) for i in 1:division-1]
    storage_index::I = 1
end

function update_storage!(b::FB, s::FBS)
    #Save important data from step
    step = b.current_step
    s.most_stressed_fiber[step] = 1/b.max_σ
    s.nr_clusters[step] = b.c # The last cluster id is also the number of clusters
    s.largest_cluster[step] = maximum(b.cluster_size)
    s.largest_perimiter[step] = maximum(b.cluster_outline_length)

    if b.spanning_cluster_id != -1 && !s.spanning_cluster_has_been_found
        s.spanning_cluster_size_storage = b.cluster_size[b.spanning_cluster_id]
        s.spanning_cluster_perimiter_storage = b.cluster_outline_length[b.spanning_cluster_id]
        s.spanning_cluster_step = step
        s.spanning_cluster_has_been_found = true
    end
end


function get_fb(L, seed; α=2.0, t=0, nr="LLS", dist="Uniform", without_storage=false)
    N=L*L
    
    @assert seed!=-1 "Seed should not be negative"
    Random.seed!(seed)

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
    fb = FB{Float64, Int64}(L=L, α=α, nr=nr, x=x)
    if without_storage
        return fb
    else
        return fb, FBS{Float64, Int64}(division=10, N=N)
    end
end

function get_fb(settings::Dict, seed, without_storage=false)
    L = settings["L"]
    α = settings["a"]
    t = settings["t"]
    nr = settings["nr"]
    dist = settings["dist"]
    return get_fb(L, seed, α=α, t=t, nr=nr, dist=dist, without_storage=without_storage)
end

# Define constants
# Since they can't be typed in julia version <1.8, we can't use this
#ALIVE = -1 #::Int64 = -1 # A fiber that has not yet been broken
#CURRENT_BORDER = -2 #::Int64 = -2 # The border of the current cluster being explored
#PAST_BORDER = -3#::Int64 = -3 # The border of a cluster that has been explored
#BROKEN = 0#::Int64 = 0 # A broken fiber that has not been explored
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

function healBundle!(b::FB; reset_break_sequence=true)
    # This completely resets the bundle
    fill!(b.σ, 1)
    fill!(b.tension, 0)
    fill!(b.status, -1)
    b.current_step = 0
    if reset_break_sequence
        fill!(b.break_sequence, 0)
    end
    fill!(b.cluster_size, 0)
    resetBundle!(b)
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
    update_σ!(b)
    for i in eachindex(b.σ)
        b.tension[i] = b.σ[i] / b.x[i]
    end
end

function find_next_fiber!(b::FB)
    b.current_step += 1
    b.break_sequence[b.current_step] = argmax(b.tension)
    b.max_σ = b.tension[b.break_sequence[b.current_step]]
end

function findAndBreakNextFiber!(b::FB, s::FBS)
    update_tension!(b)
    find_next_fiber!(b)
    update_storage!(b, s)
    break_fiber!(b)
    resetBundle!(b)
end

function findAndBreakNextFiber!(b::FB)
    findAndBreakNextFiber!(b, FBS{Float64, Int64}(division=10, N=b.N))
end

function break_fiber!(b::FB)
    i = b.break_sequence[b.current_step]
    break_this_fiber!(i, b)
end

function break_this_fiber!(i::Int64, b::FB)
    b.status[i] = 0#BROKEN
    b.σ[i] = 0
    b.tension[i] = 0
end

function break_fiber_list!(I::AbstractArray{Int}, b::FB)
    for i in I
        if i==0 # We skipp 0 for convenience.
            continue
        end
        break_this_fiber!(i, b)
    end
end


function reset_relative_possition!(b::FB)
    # Note about the relative possition:
    # If a fiber boarders two clusters, it will have the relative
    # possition of the cluster explored last.
    fill!(b.rel_pos_x, 0)
    fill!(b.rel_pos_y, 0)
end

function update_σ!(b::FB)
    @assert b.c == 0 "The bundle have not been reset!"
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

function fiber_index_to_xy(i::Int64, L::Int64)
    y = mod1(i,L)
    x = ceil(Int64, i/L) #ceil to use one indexing
    return x, y
end

function distance(a,b,L)
    # When calculating various propperties of a cluster
    # we often want the distance between two points. 
    # because of periodic boundry conditions, this
    # is slightly more difficult than in normal situations

    # This function finds the shortest distance between two NUMBERS. NOT POINTS!
    # i.e. It is used on a 1D line
    # Ex: distance (1,9,10) is not 8 like in the normal case, but 
    # 2. Because of Julias one indexing, we use mod1, and the distance
    # is therefore not 3.

    return minimum([abs(a-b), abs(a-(b-L)), abs(a-(b+L))])
end

function save_initial_cluster_possition(i::Int64, b::FB)
    # Calculating the cm is a bit difficult because
    # of the periodic boundryconditions. In order to save
    # some variables, we are going to store the initial 
    # entry possition of the cluser in the last entry of b.cluster_cm_x/y.
    # Since we can never have N clusters, this space is never
    # used and we can freely use it for what we want. 
    # Then we will use relative possitions with respect to this
    # value in further calculations instead of potentially
    # looping over to a periodic site.
    x,y = fiber_index_to_xy(i, b.L)
    b.cluster_cm_x[b.N] = x
    b.cluster_cm_y[b.N] = y
    b.rel_pos_x[i] = 0
    b.rel_pos_y[i] = 0
end

function add_to_cm(i::Int64, b::FB)
    # First get source possition
    x = b.cluster_cm_x[b.N]
    y = b.cluster_cm_y[b.N]
    # Now add that together with the relative possition
    b.cluster_cm_x[b.c] += x + b.rel_pos_x[i]
    b.cluster_cm_y[b.c] += y + b.rel_pos_y[i]
end

function normalize_cm(b::FB)
    b.cluster_cm_x[b.c] /= b.cluster_size[b.c]
    b.cluster_cm_y[b.c] /= b.cluster_size[b.c]

    # We don't want the center of mass to end up outside of the grid
    b.cluster_cm_x[b.c] = mod1(b.cluster_cm_x[b.c], b.L)
    b.cluster_cm_y[b.c] = mod1(b.cluster_cm_y[b.c], b.L)

    # Also, to not get confused, we reset
    # the initial possition so that we don't think that it is a cm
    b.cluster_cm_x[b.N]=0
    b.cluster_cm_y[b.N]=0
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
    # We set the center of mass to be zero
    b.cluster_cm_x[b.c] = 0
    b.cluster_cm_y[b.c] = 0
    # We save the initial possition
    save_initial_cluster_possition(i, b)
    # We reset cluster dimisions
    reset_cluster_dimensions!(b)

    # While there are still unexplored fibers in the cluster
    while nr_unexplored > nr_explored

        #@logmsg clusterLog "Explored $nr_unexplored / $nr_explored"
        # Preemptively count this fiber as explored (because 1 indexing)
        nr_explored += 1
        # Get the next unexplored fiber
        current_fiber = b.unexplored[nr_explored]
        # Add to the center of mass of the cluster
        add_to_cm(current_fiber, b)
        # Go through all neighbours of the fiber
        nr_unexplored = check_neighbours!(current_fiber, nr_unexplored, b)
    end
    
    # Now that we know how large the cluster is
    # we can normalize the size of the cm
    normalize_cm(b)

    return spanning(b)


end

function spanning(b::FB)
    # Cluster dimensions
    # [max_x, min_x, max_y, min_y]
    # L-1 because the relative coordinates in the cluster start at 0,0
    # NB! Once the cluster is spanning, the dimension is no longer reliable because of
    # periodicity.
    # NNB!! We actually do L-2 to prevent the cluster from connecting with itself.
    # This creates multiple center of masses! With different radii of gyration!!
    # Data before jan. 4 2023 will not have spanning clusters with size L-2.
    return abs(b.cluster_dimensions[1] - b.cluster_dimensions[2]) >= b.L-2 ||
           abs(b.cluster_dimensions[3] - b.cluster_dimensions[4]) >= b.L-2
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
    # Run over cluster outline fibers
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
    if b.nr == "LLS"
        # With the Uniform neighbourhood rule, we can apply a simple stress
        apply_simple_stress!(b)
    elseif b.nr == "ELS"
        apply_ELS_stress!(b)
    else
        # But with more complex rules, we need to do it in two steps
        # First a calculation to find the fiber strengths (As a function of their neighbourhood), and then apply the stress
        if b.nr == "CNR"
            apply_to_neighbourhood!(neighbourhoodToStrength, b)
        elseif b.nr == "CLS"
            apply_to_neighbourhood!(alive_fibers_in_neighbourhood, b)
        else
            #@debug "Unknown neighbourhood rule: $neighbourhood_rule"
            error("Unknown neighbourhood rule")
        end

        apply_stress!(b)
    end
end

function apply_ELS_stress!(b::FB)
    stress = b.N/(b.N-b.current_step)
    if b.c == 1
        for i in eachindex(b.status)
            if b.status[i] < 0 #Not Broken
                b.σ[i] = stress
            end
        end
    else
        # We only add stress to the whole bundle once
    end
end

function apply_simple_stress!(b::FB)
    added_stress = b.cluster_size[b.c]/b.cluster_outline_length[b.c]
    for i in 1:b.cluster_outline_length[b.c]
        fiber = b.cluster_outline[i]
        b.σ[fiber] += added_stress
        b.status[fiber] = -3 #PAST_BORDER
    end
end


function apply_stress!(b::FB)
    # See page 26 in Jonas Tøgersen Kjellstadli's doctoral theses, 2019:368
    # High alpha means that having neighbours is more important
    C=0.0
    for i in 1:b.cluster_outline_length[b.c]
        C += b.neighbourhood_values[i] ^(-b.α+1)
    end
    C = 1/C # A normalization constant
    for i in 1:b.cluster_outline_length[b.c]
        fiber = b.cluster_outline[i]
        g = C * b.neighbourhood_values[i] ^(-b.α+1)
        added_stress =  b.cluster_size[b.c]*g
        b.σ[fiber] += added_stress
        b.status[fiber] = -3 #PAST_BORDER
    end
end


#@time fb, storage = get_fb(100)
#pprof(;webport=58699)