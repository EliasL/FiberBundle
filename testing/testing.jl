using ProgressMeter
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
    status = fill(-1, N)
    cluster_size = zeros(Int64, N)
    cluster_dimensions = zeros(Int64, N)
    cluster_outline = zeros(Int64, N)
    cluster_outline_length = zeros(Int64, N)
    unexplored = zeros(Int64, N)

    # First fiber!
    i = findNextFiber(σ, x)
    @assert i==1 "The first fiber should break"
    resetClusters(status, σ)
    break_fiber(i, status, σ)
    @assert status[i]==BROKEN "The first fiber should be broken"
    update_σ(status, σ, neighbours, neighbourhoods, cluster_size, cluster_dimensions, cluster_outline, cluster_outline_length, unexplored)
    @assert neighbours[i, :] == [3,2,4,7] "Incorrect neighbours"
    @assert neighbours[9, :] == [8,7,3,6] "Incorrect neighbours"
    @assert all(status[neighbours[1,:]] .== PAST_BORDER) "These should be past borders"
    @assert all(σ[neighbours[1,:]] .== 1.25) "The tension is incorrect"
    @assert sum(σ) ≈ N "No conservation of tension"

    # Second fiber
    i = findNextFiber(σ, x)
    @assert i==2 "The second fiber should break"
    resetClusters(status, σ)
    break_fiber(i, status, σ)
    update_σ(status, σ, neighbours, neighbourhoods, cluster_size, cluster_dimensions, cluster_outline, cluster_outline_length, unexplored)
    @assert all(status[[1,2]] .== 1) "Both should belong to the same cluster"
    @assert sum(σ) ≈ N "No conservation of tension"

    # Break the rest for good measure
    for _ in 3:N-1
        i = findNextFiber(σ, x)
        resetClusters(status, σ)
        break_fiber(i, status, σ)
        update_σ(status, σ, neighbours, neighbourhoods, cluster_size, cluster_dimensions, cluster_outline, cluster_outline_length, unexplored)
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
        status = fill(-1, N)
        cluster_size = zeros(Int64, N)
        cluster_dimensions = zeros(Int64, N)
        cluster_outline = zeros(Int64, N)
        cluster_outline_length = zeros(Int64, N)
        unexplored = zeros(Int64, N)

        for i in broken
            break_fiber(i, status, σ)
        end
        update_σ(status, σ, neighbours, neighbourhoods, cluster_size, cluster_dimensions, cluster_outline, cluster_outline_length, unexplored)
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
        σ  = ones(Float64, N) # Relative tension
        for i in broken
            break_fiber(i, status, σ)
        end

        ids = get_id_of_neighbourhoods_of_outline(status, outline, neighbourhoods)
        @assert sort(ids) == sort(correct_ids) "Unexpected neighbourhood id\nFound    $ids\nExpected $correct_ids"
    end
end

function neighbourhood_strength_test(nr)

    for _ in 1:5
        L = 4
        N = L*L # Number of fibers
        x = ones(N) # Max extension 
        σ  = ones(Float64, N) # Relative tension
        neighbours = fillAdjacent(L, NEIGHBOURS)
        neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
        status = fill(-1, N)
        cluster_size = zeros(Int64, N)
        cluster_dimensions = zeros(Int64, N)
        cluster_outline = zeros(Int64, N)
        cluster_outline_length = zeros(Int64, N)
        unexplored = zeros(Int64, N)

        for _ in 1:N-1
            i = findNextFiber(σ, x)
            resetClusters(status, σ)
            break_fiber(i, status, σ)
            update_σ(status, σ, neighbours, neighbourhoods, cluster_size, cluster_dimensions, cluster_outline, cluster_outline_length, unexplored; neighbourhood_rule=nr)
            @assert sum(σ) ≈ N "No conservation of tension"
        end
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
        neighbourhood_strength_test(nr)
    end
    println("Neighbourhood strength test complete")
    println("All tests completed!")
end

test()