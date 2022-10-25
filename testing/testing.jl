using ProgressMeter
using Plots

include("../burningMan.jl")



seed=0
Random.seed!(seed) # Setting the seed

function basic_test()
    # This test will break the fibers in this order
    # 1 4 7
    # 2 5 8
    # 3 6 9

    L = 3
    N = L*L # Number of fibers
    x = (1:N)./N # Max extension 
    σ  = ones(Float64, N) # Relative tension
    neighbours = fillAdjacent(L, NEIGHBOURS)
    neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
    neighbourhood_values = zeros(Int64, N)
    status = fill(-1, N)
    cluster_size = zeros(Int64, N)
    cluster_dimensions = zeros(Int64, 4)
    rel_pos_x = zeros(Int64, N)
    rel_pos_y = zeros(Int64, N)
    cluster_outline = zeros(Int64, N)
    cluster_outline_length = zeros(Int64, N)
    unexplored = zeros(Int64, N)

    # First fiber!
    i = findNextFiber(σ, x)
    @assert i==1 "The first fiber should break"
    resetClusters(status, σ)
    break_fiber(i, status, σ)
    @assert status[i]==BROKEN "The first fiber should be broken"
    update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored)
    @assert neighbours[i, :] == [4,7,3,2] "Incorrect neighbours: $(neighbours[i, :])"
    @assert neighbours[9, :] == [3,6,8,7] "Incorrect neighbours: $(neighbours[9, :])"
    @assert all(status[neighbours[1,:]] .== PAST_BORDER) "These should be past borders"
    @assert all(σ[neighbours[1,:]] .== 1.25) "The tension is incorrect"
    @assert sum(σ) ≈ N "No conservation of tension"

    # Second fiber
    i = findNextFiber(σ, x)
    @assert i==2 "The second fiber should break"
    resetClusters(status, σ)
    break_fiber(i, status, σ)
    update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored)
    @assert all(status[[1,2]] .== 1) "Both should belong to the same cluster"
    @assert sum(σ) ≈ N "No conservation of tension"

    # Break the rest for good measure
    for _ in 3:N-1
        i = findNextFiber(σ, x)
        resetClusters(status, σ)
        break_fiber(i, status, σ)
        update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored)
        @assert sum(σ) ≈ N "No conservation of tension"
    end
end

function cluster_test()
    #  1  5  9 13
    #  2  6 10 14
    #  3  7 11 15
    #  4  8 12 16
    tests = [
        [
            [6,7,10,11], #broken
            4, # Cluster size
            [2,3,5,8,9,12,14,15], # Outline
            8 # Outline length
        ],
        [
            [1,2,5], #broken
            3, # Cluster size
            [3,4,6,9,8,13,14], # Outline
            7 # Outline length
        ],
    ]

    for sol in tests
        broken, clusterSize, outline, length = sol
        L = 4
        N = L*L # Number of fibers
        x = ones(N) # Max extension 
        σ  = ones(Float64, N) # Relative tension
        neighbours = fillAdjacent(L, NEIGHBOURS)
        neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
        neighbourhood_values = zeros(Int64, N)
        status = fill(-1, N)
        cluster_size = zeros(Int64, N)
        cluster_dimensions = zeros(Int64, 4)
        rel_pos_x = zeros(Int64, N)
        rel_pos_y = zeros(Int64, N)
        cluster_outline = zeros(Int64, N)
        cluster_outline_length = zeros(Int64, N)
        unexplored = zeros(Int64, N)

        for i in broken
            break_fiber(i, status, σ)
        end
        update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored)
        @assert sum(σ) ≈ N "No conservation of tension"

        @assert clusterSize == cluster_size[1] "Cluster size missmatch"
        @assert length == cluster_outline_length[1] "Cluster outline length missmatch"
        @assert sort(outline) == sort(cluster_outline[1:length]) "Cluster outline missmatch"
    end
end

function neighbourhood_id_test()
    #  1  5  9 13
    #  2  6 10 14
    #  3  7 11 15
    #  4  8 12 16
    tests = [
        [
            4, #L
            [7], #broken
            [3,6,8,11], # Outline
            [240, 192, 254, 248], # Correct ids
        ],
    ]
    for sol in tests
        L, broken, outline, correct_ids = sol
        N = L*L
        status = fill(-1, N)
        neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
        neighbourhood_values = zeros(Int64, N)
        σ  = ones(Float64, N) # Relative tension
        for i in broken
            break_fiber(i, status, σ)
        end

        apply_to_neighbourhood(neighbourhoodToInt, status, outline, length(outline), neighbourhood_values, neighbourhoods)
        @assert sort(neighbourhood_values[1:length(outline)]) == sort(correct_ids) "Unexpected neighbourhood id\nFound    $neighbourhood_values\nExpected $correct_ids"
    end
end

function neighbourhood_strength_test(nr)

    for _ in 1:5
        L = 4
        N = L*L # Number of fibers
        x = rand(N) # Max extension 
        σ  = ones(Float64, N) # Relative tension
        neighbours = fillAdjacent(L, NEIGHBOURS)
        neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
        neighbourhood_values = zeros(Int64, N)
        status = fill(-1, N)
        cluster_size = zeros(Int64, N)
        cluster_dimensions = zeros(Int64, 4)
        rel_pos_x = zeros(Int64, N)
        rel_pos_y = zeros(Int64, N)
        cluster_outline = zeros(Int64, N)
        cluster_outline_length = zeros(Int64, N)
        unexplored = zeros(Int64, N)

        for _ in 1:N-1
            i = findNextFiber(σ, x)
            resetClusters(status, σ)
            break_fiber(i, status, σ)
            update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size,
                     cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length,
                     unexplored; neighbourhood_rule=nr,α=2.0)
            @assert sum(σ) ≈ N "No conservation of tension.\n$(sum(σ)) ≠ $N "
        end
    end
end

function neighbourhood_strength_test_with_alpha(nr)

    for _ in 1:5
        L = 4
        N = L*L # Number of fibers
        x = rand(N) # Max extension 
        α = 1+rand()*10
        σ  = ones(Float64, N) # Relative tension
        neighbours = fillAdjacent(L, NEIGHBOURS)
        neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
        neighbourhood_values = zeros(Int64, N)
        status = fill(-1, N)
        cluster_size = zeros(Int64, N)
        cluster_dimensions = zeros(Int64, 4)
        rel_pos_x = zeros(Int64, N)
        rel_pos_y = zeros(Int64, N)
        cluster_outline = zeros(Int64, N)
        cluster_outline_length = zeros(Int64, N)
        unexplored = zeros(Int64, N)

        for _ in 1:N-1
            i = findNextFiber(σ, x)
            resetClusters(status, σ)
            break_fiber(i, status, σ)
            
            update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size,
                     cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length,
                     unexplored; neighbourhood_rule=nr,α=α)
            @assert sum(σ) ≈ N "No conservation of tension.\n$(sum(σ)) ≠ $N "
        end
    end
end

function test_store_possition()
    current_fiber = 1
    neighbour_fiber = 4
    # Directions = [+1x,-1x,+1y,-1y]
    direction = 1 #x+1
    cluster_dimensions= [2,-4,1,1]
    rel_pos_x = [-4,2,0,0]
    rel_pos_y= [1,1,0,0]

    store_possition(current_fiber, neighbour_fiber, direction, cluster_dimensions, rel_pos_x, rel_pos_y)
    @assert rel_pos_x[neighbour_fiber] == -3 "x movement not working"
    @assert cluster_dimensions[2] == -4 "Cluster dimensions not working1: $cluster_dimensions"

    direction = 2 #x-1
    store_possition(current_fiber, neighbour_fiber, direction, cluster_dimensions, rel_pos_x, rel_pos_y)
    @assert rel_pos_x[neighbour_fiber] == -5 "x movement not working"
    @assert cluster_dimensions[2] == -5 "Cluster dimensions not working1: $cluster_dimensions"

    direction = 3 #y+1
    store_possition(current_fiber, neighbour_fiber, direction, cluster_dimensions, rel_pos_x, rel_pos_y)
    @assert cluster_dimensions == [2,-5,2,1] "Cluster dimensions not working2: $cluster_dimensions"
    @assert rel_pos_y[neighbour_fiber] == 2 "y movement not working1: $rel_pos_y"

    direction = 4 #y-1
    store_possition(current_fiber, neighbour_fiber, direction, cluster_dimensions, rel_pos_x, rel_pos_y)
    @assert cluster_dimensions == [2,-5,2,0] "Cluster dimensions not working3: $cluster_dimensions"
    @assert rel_pos_y[neighbour_fiber] == 0 "y movement not working2: $rel_pos_y"
end


function spanning_cluster_test()

    L = 16
    N = L*L # Number of fibers
    x = rand(N) # Max extension 
    σ  = ones(Float64, N) # Relative tension
    neighbours = fillAdjacent(L, NEIGHBOURS)
    neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
    neighbourhood_values = zeros(Int64, N)
    status = [0, -1, -1, -1, 0, -1, -1, -1, 0, -1, 0, 0, -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, 0, 0, 0, 0, -1, -1, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0]
    
    cluster_size = zeros(Int64, N)
    cluster_dimensions = zeros(Int64, 4)
    rel_pos_x = zeros(Int64, N)
    rel_pos_y = zeros(Int64, N)
    cluster_outline = zeros(Int64, N)
    cluster_outline_length = zeros(Int64, N)
    unexplored = zeros(Int64, N)
    spanning_cluster_size = 0
    c, spanning_cluster, spanning_cluster_size = update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored)

    if spanning_cluster != -1
        @assert spanning_cluster_size >= L "Impossibly small cluster"
        
    end 
end

function random_spanning_cluster_test()

    for run_nr in 1:5
        L = 32
        N = L*L # Number of fibers
        x = rand(N) # Max extension 
        σ  = ones(Float64, N) # Relative tension
        neighbours = fillAdjacent(L, NEIGHBOURS)
        neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
        neighbourhood_values = zeros(Int64, N)
        status = fill(-1, N)
        cluster_size = zeros(Int64, N)
        cluster_dimensions = zeros(Int64, 4)
        rel_pos_x = zeros(Int64, N)
        rel_pos_y = zeros(Int64, N)
        cluster_outline = zeros(Int64, N)
        cluster_outline_length = zeros(Int64, N)
        unexplored = zeros(Int64, N)
        spanning_cluster_size = 0

        for k in 1:N
            i = findNextFiber(σ, x)
            resetClusters(status, σ)
            break_fiber(i, status, σ)
            c, spanning_cluster, spanning_cluster_size = update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored)

            if spanning_cluster != -1
                break
                @assert spanning_cluster_size >= L "Impossibly small cluster"
                
            end 
        end
        #@assert cluster_dimensions[1] - cluster_dimensions[2] == cluster_dimensions[3] - cluster_dimensions[4] == L-1 "Not correct dimensions! $cluster_dimensions"
        #@assert spanning_cluster_size == N "Not correct cluster size:  $spanning_cluster_size ≠ $N "
    end
end

function test()
    basic_test()
    println("Basic test complete")
    cluster_test()
    println("Cluster test complete")
    neighbourhood_id_test()
    println("Neighbourhood id test complete")
    for nr in ["UNR", "CNR", "SNR"]
        print(nr*"\r")
        neighbourhood_strength_test(nr)
    end
    println("Neighbourhood strength test complete")
    test_store_possition()
    spanning_cluster_test()
    random_spanning_cluster_test()
    println("Cluster dimensions tests complete")
    println("All tests completed!")
end

test()
#spanning_cluster_test()
#random_spanning_cluster_test()