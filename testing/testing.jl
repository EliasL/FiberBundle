using ProgressMeter
using Plots
using Test

include("../burningMan.jl")
include("../plottingScripts/showBundle.jl")
include("../dataGenerator.jl")
include("../support/inertia.jl")

seed=0
Random.seed!(seed) # Setting the seed

logger = SimpleLogger(stdout, 10000)
global_logger(logger)


function basic_test() 
    # This test will break the fibers in this order
    # 1 4 7
    # 2 5 8
    # 3 6 9
    L=3
    b = get_fb(L, dist= f(N)=(1:N)./N, without_storage=true) #set distribution to break in fixed order
    # First fiber!
    @test sum(b.status) == -b.N #All fibers are alive
    findNextFiber!(b)
    i = b.break_sequence[b.current_step]
    @test i==1 #"The first fiber should break"
    break_fiber!(b)
    @test b.status[i]== 0#BROKEN #"The first fiber should be broken"
    @test sum(b.status) == (-b.N +1) # Now all but one fiber should be broken
    update_σ!(b)
    @test b.neighbours[i, :] == [4,7,3,2] #"Incorrect neighbours: $(neighbours[i, :])"
    @test b.neighbours[9, :] == [3,6,8,7] #"Incorrect neighbours: $(neighbours[9, :])"
    @test all(b.status[b.neighbours[1,:]] .== -3)#PAST_BORDER) #"These should be past borders"
    @test all(b.σ[b.neighbours[1,:]] .== 1.25) #"The tension is incorrect"
    @test sum(b.σ) ≈ b.N #"No conservation of tension"
    
    # Second fiber
    findNextFiber!(b)
    @test b.break_sequence[b.current_step]==2 #"The second fiber should break"
    resetBundle!(b)
    break_fiber!(b)
    update_σ!(b)
    @test all(b.status[[1,2]] .== 1) #"Both should belong to the same cluster"
    @test sum(b.σ) ≈ b.N #"No conservation of tension"
    
    # Break the rest for good measure
    for _ in 3:b.N-1
        findNextFiber!(b)
        resetBundle!(b)
        break_fiber!(b)
        update_σ!(b)
        @test sum(b.σ) ≈ b.N #"No conservation of tension"
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
    for α in [1.0, 2, 2.2, 3.3]
        L = 4
        b = get_fb(L, α=α, nr=nr, without_storage=true)

        for _ in 1:b.N-1
            findNextFiber!(b)
            resetBundle!(b)
            break_fiber!(b)
            update_σ!(b)
            @test sum(b.σ) ≈ b.N #"No conservation of tension.\n$(sum(σ)) ≠ $N "
        end
    end
end

function test_store_possition()
    current_fiber = 1
    neighbour_fiber = 4
    # Directions = [+1x,-1x,-1y,1y]
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
    @test b.cluster_dimensions == [2,-5,0,1] #"Cluster dimensions not working2: $cluster_dimensions"
    @test b.rel_pos_y[neighbour_fiber] == 0 #"y movement not working1: $rel_pos_y"

    direction = 4 #y-1
    store_possition!(current_fiber, neighbour_fiber, direction, b)
    @test b.cluster_dimensions == [2,-5,0,2] #"Cluster dimensions not working3: $cluster_dimensions"
    @test b.rel_pos_y[neighbour_fiber] == 2 #"y movement not working2: $rel_pos_y"
end


function spanning_cluster_test()
    L = 16
    b = get_fb(L, without_storage=true)
    test_s = [0, -1, -1, -1, 0, -1, -1, -1, 0, -1, 0, 0, -1, 0, -1, -1, 0, -1, -1, 0, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, 0, 0, 0, 0, -1, -1, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, 0, 0, -1, -1, 0, 0, 0, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, 0, 0, -1, -1, -1, 0, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1, -1, 0, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, 0, -1, -1, -1, -1, 0, -1, -1, -1, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0, -1, -1, -1, -1, -1, -1, 0, -1, 0, -1, 0, 0, 0, 0]
    b.status .= test_s
    @test b.spanning_cluster_id == -1 #Not spanning
    update_σ!(b)
    @test b.spanning_cluster_id != -1 #Spanning!
    @test b.cluster_size[b.spanning_cluster_id] == 74 #This is probably correct
    #plot_fb(b)
    
    L = 16
    b = get_fb(L, without_storage=true)
    m = reshape(test_s, (L,L))
    m = rotl90(m)
    b.status = reshape(m, (L*L))
    @test b.spanning_cluster_id == -1 #Not spanning
    update_σ!(b)
    #plot_fb(b)
    @test b.spanning_cluster_id != -1 #Spanning!
    @test b.cluster_size[b.spanning_cluster_id] == 74 #This is probably correct    
    
    L = 16
    b = get_fb(L, without_storage=true)
    m = reshape(test_s, (L,L))
    m = rot180(m)
    b.status = reshape(m, (L*L))
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
            findNextFiber!(b)
            resetBundle!(b)
            break_fiber!(b)
            update_σ!(b)

            if b.spanning_cluster_id != -1
                @test b.cluster_size[b.spanning_cluster_id] >= L #"Impossibly small cluster"
                break
            end 
        end
    end
end

function basic_cm_test()
    # 1,4,7
    # 2,5,8
    # 3,6,9
    L=3
    b = get_fb(L, without_storage=true)
    for i in [1,2,4,5]
        break_this_fiber!(i,b)
    end
    update_σ!(b)
    #@test b.cluster_cm_x[9] == 1
    #@test b.cluster_cm_y[9] == 1
    @test b.cluster_cm_x[1] == 1.5
    @test b.cluster_cm_y[1] == 1.5


    b = get_fb(L, without_storage=true)
    for i in [4,5,6,7,8,9]
        break_this_fiber!(i,b)
    end
    update_σ!(b)
    #@test b.cluster_cm_x[9] == 2
    #@test b.cluster_cm_y[9] == 1
    @test b.cluster_cm_x[1] == 2.5
    # This is a very important test. The cm of this cluster
    # can actually have a y value of 1, 2 or 3. Since the boundry
    # is periodic, it is onlt a question of reference
    @test b.cluster_cm_y[1] in [1,2,3]
end

function periodic_cm_test()
    
    #  1  5  9 13
    #  2  6 10 14
    #  3  7 11 15
    #  4  8 12 16
    L=4
    b = get_fb(L, without_storage=true)
    b.status = [-1,  -1,  0,  -1,
    -1,  -1, -1,  -1,
    -1,  -1,  0,  -1,
    -1,  -1,  0,  -1]
    update_σ!(b)
    @test b.cluster_cm_x[1] == 4
    @test b.cluster_cm_y[1] == 3
    
    #plot_fb(b, show=false, axes=true)
    #display(plot_fb_cm(b))
end

function periodic_distance_test()
    @test distance(0,0,1) == 0
    @test distance(1,0,5) == 1
    @test distance(0,1,5) == 1
    @test distance(1,4,5) == 2
    @test distance(4,1,5) == 2
    @test distance(4,5,5) == 1
    @test distance(5,4,5) == 1
    @test distance(8,4,5) == 1
    @test distance(4,8,5) == 1
end

function storageTest()
    test_data_path = "test_data/"
    settings = make_settings(8, 0.1, "CLS", 2.0, test_data_path)

    # Test clean generation
    seeds = 1:3
    for seed in seeds
        break_bundle(settings, nothing, nothing, seed, use_threads=false, stop_after_spanning=false, use_past_progress=false)
        @test ispath(get_file_name(settings, seed, false))
    end
    clean_after_run(settings, seeds)
    @test ispath(get_file_name(settings, -1, false))
    @test ispath(get_file_name(settings, -1, true))

    # Test overwrite and only go to spanning
    break_bundle(settings, nothing, nothing, 1, use_threads=false)
    clean_after_run(settings, [1])
    @test ispath(get_file_name(settings, -1, true))

    seed = 1
    f = load_file(settings, seed=-1, average=false)
    for key in data_keys
        @test key*"/$seed" in keys(f)
    end
    f = load_file(settings, seed=-1, average=true)
    for key in averaged_data_keys
        @test "average_"*key in keys(f)
    end
    expand_file(settings)
    f = load_file(settings, seed=3, average=false)
    for key in data_keys
        @test key in keys(f)
    end


    # Test to see if a bundle can partially be broken and continue later

    settings = make_settings(32, 0.1, "CLS", 2.0, test_data_path)
    # First do everything at once to create a correct answer
    break_bundle(settings, nothing, nothing, 1, use_threads=false,stop_after_spanning=false, use_past_progress=false)
    clean_after_run(settings, [1])
    correct_f = load_file(settings, average=false)

    # Now do it splitt and see if we can reproduce the same result
    break_bundle(settings, nothing, nothing, 1, use_threads=false,stop_after_spanning=true, use_past_progress=false)
    break_bundle(settings, nothing, nothing, 1, use_threads=false,stop_after_spanning=false, use_past_progress=true)
    clean_after_run(settings, [1])
    test_f = load_file(settings, average=false)

    @assert all(correct_f["largest_cluster/$seed"] .== test_f["largest_cluster/$seed"])
    @assert abs(correct_f["simulation_time/$seed"] - test_f["simulation_time/$seed"]) < correct_f["simulation_time/$seed"] *0.05



    search_for_loose_files(settings)
    rm(test_data_path, force=true, recursive=true)
end

function find_strange_fb()
    for run_nr in 1:500
        L = 4
        Random.seed!(run_nr)
    
        b=get_fb(L, without_storage=true)

        for k in 1:b.N
            findNextFiber!(b)
            resetBundle!(b)
            break_fiber!(b)
            update_σ!(b)
            if !all(0 .<= b.cluster_cm_x .<= b.L) || !all(0 .<= b.cluster_cm_y .<= b.L)
                println(b.cluster_cm_x)
                println(b.cluster_cm_y)
                resetBundle!(b)
                println(b.status)
                plot_fb(b, show=false, axes=true)
                display(plot_fb_cm(b))
                return 
            end
        end
    end
end
#find_strange_fb()


function intertiaTest()
    #  1  5  9 13
    #  2  6 10 14
    #  3  7 11 15
    #  4  8 12 16
    L=4
    b = get_fb(L, without_storage=true)
    b.status = [-1,  -1,  0,  -1,
                -1,  -1, -1,  -1,
                -1,  -0,  0,  -1,
                -1,  -1,  0,  -1]
    update_σ!(b)

    minor_axes, major_axes, minor_values, major_values = find_major_and_minor_axes(b)
    p = plot_fb(b, show=false)
    plot_fb_axes(b, minor_axes, major_axes, minor_values, major_values)
    display(p)
end

#intertiaTest()

function custom_cluster_test()
    #  1  5  921 13
    #  2  6 10 14
    #  3  7 11 15
    #  4  8 12 16
    broken = [6,10,11]
        
    L=4
    b = get_fb(L, nr="CLS", α=0.5, dist=ones(L*L), without_storage=true) #set distribution to break in fixed order

    for i in broken
        break_this_fiber!(i, b)
    end
    update_σ!(b)
    display(reshape(b.σ, (L,L)))
end
#custom_cluster_test()

function test()
    
    @testset verbose=true "Tests" begin
        
        @testset "Basic" begin basic_test() end
        @testset "Cluster" begin cluster_test() end
        @testset "Neighbourhood id" begin neighbourhood_id_test() end
        @testset "Neighbourhood rules $nr" for nr in ["LLS", "CNR", "CLS"]
            neighbourhood_strength_test_with_alpha(nr)
        end
        @testset "Store possition" begin test_store_possition() end
        @testset "Spanning cluster" begin spanning_cluster_test() end
        @testset "Ransom spanning cluster" begin random_spanning_cluster_test() end
        @testset "Center of mass" begin basic_cm_test() 
                                        periodic_cm_test() end
        @testset "Distance" begin periodic_distance_test() end
        #@testset "1D cluster size" begin end

        @testset "Storage" begin storageTest() end
        
    end
end
test()
println("Test done")