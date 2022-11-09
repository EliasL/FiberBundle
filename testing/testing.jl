using ProgressMeter
using Plots
using Test

include("../burningMan.jl")


seed=0
Random.seed!(seed) # Setting the seed

function basic_test()
    @testset "Basic tests" begin       
        # This test will break the fibers in this order
        # 1 4 7
        # 2 5 8
        # 3 6 9
        L=3
        b = get_fb(L, dist= f(N)=(1:N)./N, without_storage=true) #set distribution to break in fixed order
        # First fiber!
        @test sum(b.status) == -b.N #All fibers are alive
        i = findNextFiber!(b)
        @test i==1 #"The first fiber should break"
        break_this_fiber!(i, b)
        @test b.status[i]==BROKEN #"The first fiber should be broken"
        @test sum(b.status) == (-b.N +1) # Now all but one fiber should be broken
        update_σ!(b)
        @test b.neighbours[i, :] == [4,7,3,2] #"Incorrect neighbours: $(neighbours[i, :])"
        @test b.neighbours[9, :] == [3,6,8,7] #"Incorrect neighbours: $(neighbours[9, :])"
        @test all(b.status[b.neighbours[1,:]] .== PAST_BORDER) #"These should be past borders"
        @test all(b.σ[b.neighbours[1,:]] .== 1.25) #"The tension is incorrect"
        @test sum(b.σ) ≈ b.N #"No conservation of tension"
        
        # Second fiber
        i = findNextFiber!(b)
        @test i==2 #"The second fiber should break"
        resetBundle!(b)
        break_this_fiber!(i, b)
        update_σ!(b)
        @test all(b.status[[1,2]] .== 1) #"Both should belong to the same cluster"
        @test sum(b.σ) ≈ b.N #"No conservation of tension"
        
        # Break the rest for good measure
        for _ in 3:b.N-1
            i = findNextFiber!(b)
            resetBundle!(b)
            break_this_fiber!(i, b)
            update_σ!(b)
            @test sum(b.σ) ≈ b.N #"No conservation of tension"
        end
    end
end

function cluster_test()
    #  1  5  921 13
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

    @testset "Cluster test #$i" for (i,sol) in enumerate(tests)
        broken, clusterSize, outline, length = sol
        
        L=4
        b = get_fb(L, dist=ones(L*L), without_storage=true) #set distribution to break in fixed order

        for i in broken
            break_this_fiber!(i, b)
        end
        update_σ!(b)
        @test sum(b.σ) ≈ b.N #"No conservation of tension"

        @test clusterSize == b.cluster_size[1] #"Cluster size missmatch"
        @test length == b.cluster_outline_length[1] #"Cluster outline length missmatch"
        @test sort(outline) == sort(b.cluster_outline[1:length]) #"Cluster outline missmatch"
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
        b = get_fb(L, without_storage=true)
        for i in broken
            break_this_fiber!(i,b)
        end
        update_σ!(b)

        apply_to_neighbourhood!(arr_to_int, b)
        @test sort(b.neighbourhood_values[1:length(outline)]) == sort(correct_ids) #"Unexpected neighbourhood id\nFound    $neighbourhood_values\nExpected $correct_ids"
    end
end

function neighbourhood_strength_test_with_alpha(nr)
    for _ in 1:5
        L = 4
        b = get_fb(L, α= 1+rand()*10, nr=nr, without_storage=true)

        for _ in 1:b.N-1
            i = findNextFiber!(b)
            resetBundle!(b)
            break_this_fiber!(i, b)
            update_σ!(b)
            @test sum(b.σ) ≈ b.N #"No conservation of tension.\n$(sum(σ)) ≠ $N "
        end
    end
end

function test_store_possition()
    current_fiber = 1
    neighbour_fiber = 4
    # Directions = [+1x,-1x,+1y,-1y]
    direction = 1 #x+1
    b = get_fb(2, without_storage=true)
    b.cluster_dimensions= [2,-4,1,1]
    b.rel_pos_x = [-4,2,0,0]
    b.rel_pos_y= [1,1,0,0]

    store_possition!(current_fiber, neighbour_fiber, direction, b)
    @test b.rel_pos_x[neighbour_fiber] == -3 #"x movement not working"
    @test b.cluster_dimensions[2] == -4 #"Cluster dimensions not working1: $cluster_dimensions"

    direction = 2 #x-1
    store_possition!(current_fiber, neighbour_fiber, direction, b)
    @test b.rel_pos_x[neighbour_fiber] == -5 #"x movement not working"
    @test b.cluster_dimensions[2] == -5 #"Cluster dimensions not working1: $cluster_dimensions"

    direction = 3 #y+1
    store_possition!(current_fiber, neighbour_fiber, direction, b)
    @test b.cluster_dimensions == [2,-5,2,1] #"Cluster dimensions not working2: $cluster_dimensions"
    @test b.rel_pos_y[neighbour_fiber] == 2 #"y movement not working1: $rel_pos_y"

    direction = 4 #y-1
    store_possition!(current_fiber, neighbour_fiber, direction, b)
    @test b.cluster_dimensions == [2,-5,2,0] #"Cluster dimensions not working3: $cluster_dimensions"
    @test b.rel_pos_y[neighbour_fiber] == 0 #"y movement not working2: $rel_pos_y"
end


function spanning_cluster_test()
    L = 16
    b = get_fb(L, without_storage=true)
    b.status = [0, -1, -1, -1, 0, -1, -1, -1, 0, -1, 0, 0, -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, 0, 0, 0, 0, -1, -1, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0]
    
    @test b.spanning_cluster_id == -1 #Not spanning
    update_σ!(b)
    @test b.spanning_cluster_id != -1 #Spanning!
    @test b.cluster_size[b.spanning_cluster_id] == 74 #This is probably correct
end

function random_spanning_cluster_test()
    for run_nr in 1:5
        L = 32
        b=get_fb(L, without_storage=true)

        for k in 1:b.N
            i = findNextFiber!(b)
            resetBundle!(b)
            break_this_fiber!(i, b)
            update_σ!(b)

            if b.spanning_cluster_id != -1
                @test b.cluster_size[b.spanning_cluster_id] >= L #"Impossibly small cluster"
                break
            end 
        end
    end
end

function cm_test()
    # 1,4,7
    # 2,5,8
    # 3,6,9
    L=3
    b = get_fb(L, without_storage=true)
    for i in [1,2,4,5]
        break_this_fiber!(i,b)
    end
    update_σ!(b)
    println(b.rel_pos_x)
    println(b.rel_pos_y)
    @test b.cluster_cm_x[9] == 1
    @test b.cluster_cm_y[9] == 1
    @test b.cluster_cm_x[1] == 1.5
    @test b.cluster_cm_y[1] == 1.5


    b = get_fb(L, without_storage=true)
    for i in [4,5,6,7,8,9]
        break_this_fiber!(i,b)
    end
    update_σ!(b)
    display(reshape(b.status, (3,3)))
    println("hello", reshape(b.status, (3,3))[2,1])
    println(b.rel_pos_x)
    println(b.rel_pos_y)
    @test b.cluster_cm_x[9] == 2
    @test b.cluster_cm_y[9] == 1
    @test b.cluster_cm_x[1] == 2.5
    @test b.cluster_cm_y[1] == 2
end

cm_test()
@testset verbose=true "Tests" begin
    
    @testset "Basic" begin basic_test() end
    @testset "Cluster" begin cluster_test() end
    @testset "Neighbourhood id" begin neighbourhood_id_test() end
    @testset "Neighbourhood rules $nr" for nr in ["UNR", "CNR", "SNR"]
        neighbourhood_strength_test_with_alpha(nr)
    end
    @testset "Store possition" begin test_store_possition() end
    @testset "Spanning cluster" begin spanning_cluster_test() end
    @testset "Ransom spanning cluster" begin random_spanning_cluster_test() end
    @testset "Center of mass" begin cm_test() end
end
print("Done")