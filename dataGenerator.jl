using Distributed
using ProgressMeter
using Dates

@everywhere begin
    using JLD2
end


print("Preparing workers... ")

@everywhere include("burningMan.jl")
@everywhere include("support/dataManager.jl")
@everywhere include("support/distributions.jl")

@everywhere function break_bundle(L, distribution::Function, progress_channel, working_channel, file_name, neighbourhood_rule; seed=0)
    put!(working_channel, true) # Indicate a process has started
    N = L*L # Number of fibers
    @assert seed != -1 ""
    Random.seed!(seed)
    x = distribution(N) # Max extension (Distribution)
    neighbours = fillAdjacent(L, NEIGHBOURS)
    neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)
    neighbourhood_values = zeros(Int64, N)
    # These values are reset for each step
    σ  = ones(Float64, N) # Relative tension
    max_σ = Float64(0)
    status = fill(-1,N)
    cluster_size = zeros(Int64, N)
    cluster_dimensions = zeros(Int64, 4) #Or 4? max_x, min_x, max_y, min_y
    rel_pos_x = zeros(Int64, N)
    rel_pos_y = zeros(Int64, N)
    cluster_outline_length = zeros(Int64, N)
    # These values are reset for each cluster
    cluster_outline = zeros(Int64, N)
    unexplored = zeros(Int64, N)

    # These arrays store one value for each step
    steps = N # just for readability
    most_stressed_fiber = zeros(Float64, steps)
    nr_clusters = zeros(Int64, steps)
    largest_cluster = zeros(Int64, steps)
    largest_perimiter = zeros(Int64, steps)

    # We want to store some samples of the processed
    # I'm thinking at 10%, 20%, ... 90% done would work
    # ie, 9 images
    division = 10
    if seed <= 10 #Only save samples from 10 first seeds
        status_storage = zeros(Int64, division-1, N)
        tension_storage = zeros(Float64, division-1, N)
    end
    spanning_cluster_storage = zeros(Int64, N)
    spanning_cluster_size_storage = 0
    spanning_cluster_perimiter_storage = 0
    spanning_cluster_has_not_been_found = true
    spanning_cluster_step = 0
    # If N=100 Steps to store is now [90, 80, ... , 10]
    steps_to_store = [round(Int64,N/division * i) for i in 1:division-1]
    storage_index = 1

    # Spanning cluster storage
    # Spanning cluster 
    # Break the bundle
    for step in 1:N

        # Simulate step
        i = findNextFiber(σ, x)
        max_σ = σ[i]/x[i]
        resetClusters(status, σ)
        break_fiber(i, status, σ)
        _nr_clusters, spanning_cluster, spanning_cluster_size = update_σ(status, σ, neighbours, neighbourhoods, neighbourhood_values, cluster_size, cluster_dimensions, rel_pos_x, rel_pos_y, cluster_outline, cluster_outline_length, unexplored; neighbourhood_rule=neighbourhood_rule)

        # Save important data from step
        most_stressed_fiber[step] = 1/max_σ
        nr_clusters[step] = _nr_clusters
        largest_cluster[step] = maximum(cluster_size)
        largest_perimiter[step] = maximum(cluster_outline_length)

        # Save step for visualization
        if seed <= 10 #Only save samples from 10 first seeds
            if step == steps_to_store[storage_index]
                status_storage[storage_index, :] = status
                tension_storage[storage_index, :] = σ ./ x
                if storage_index < length(steps_to_store)
                    storage_index += 1  
                end
            end
        end
        if spanning_cluster != -1 && spanning_cluster_has_not_been_found
            spanning_cluster_storage = status
            spanning_cluster_size_storage = spanning_cluster_size
            spanning_cluster_perimiter_storage = cluster_outline_length[spanning_cluster]
            spanning_cluster_step = step
            spanning_cluster_has_not_been_found = false
        end
        put!(progress_channel, true) # trigger a progress bar update
    end 
    
    @assert spanning_cluster_has_not_been_found == false "There should be a spanning cluster"

    jldopen(file_name, "w") do file
        if seed <= 10
            file["sample_states"] = status_storage
            file["tension"] = tension_storage
            file["spanning_cluster"] = spanning_cluster_storage
        end
        file["spanning_cluster_size"] = spanning_cluster_size_storage
        file["spanning_cluster_perimiter"] = spanning_cluster_perimiter_storage
        file["spanning_cluster_step"] = spanning_cluster_step
        file["sample_states_steps"] = steps_to_store./N
        file["most_stressed_fiber"] = most_stressed_fiber
        file["nr_clusters"] = nr_clusters./N
        file["largest_cluster"] = largest_cluster./N
        file["largest_perimiter"] = largest_perimiter./N
    end
    put!(working_channel, false) # trigger a progress bar update

end

# Done preparing workers!
println("Done!")


function run_workers(L, distribution_name, distribution_function, seeds, path, neighbourhood_rule)
    p = Progress(length(seeds)*L^2)
    progress = RemoteChannel(()->Channel{Bool}(), 1)
    working = RemoteChannel(()->Channel{Bool}(), 1)
    active_workers = Threads.Atomic{Int}(0)
    completed_runs = Threads.Atomic{Int}(0)


    @sync begin # start two tasks which will be synced in the very end
        # the first task updates the progress bar
        @async while take!(progress)
            ProgressMeter.next!(p; showvalues = [("Active workers", active_workers[]), ("Completed tasks", completed_runs[])])
        end

        @async while true
            started = take!(working)
            if started
                active_workers[] += 1
            else 
                active_workers[] -=1 
                completed_runs[] += 1
                if length(seeds) == completed_runs[]
                    put!(progress, false) # this tells the printing task to finish
                    break
                end
            end
        end

        # Test
        #name = get_name(L, distribution_name, path, 1)
        #break_bundle(L, distribution_function, progress, working, name, neighbourhood_rule; seed=1)

        # the second task does the computation
        @async begin
            @distributed (+) for i in seeds
                name = get_name(L, distribution_name, path, i)
                break_bundle(L, distribution_function, progress, working, name, neighbourhood_rule; seed=i)
                i^2
            end
        end
    end
end

function generate_data(path, L, requested_seeds, distribution_name, t₀, overwrite, neighbourhood_rule)

    if !isdir(path)
        println("Creating folder... ")
        mkdir(path)
    end

    missing_seeds = prepare_run(L, distribution_name, path, requested_seeds, overwrite)

    # get distribtion function
    distribution_function(n) = zeros(n) # Default

    if  occursin("Uniform", distribution_name)
        distribution_function = get_uniform_distribution(t₀)
    else
        error("No distribution found!")
    end

    if length(missing_seeds) > 0
        println("Running workers...                            ")
        run_workers(L, distribution_name, distribution_function, missing_seeds, path, neighbourhood_rule)
        println("Done!")


        print("Calculating averages... ")
        clean_after_run(L, distribution_name, path, requested_seeds)
        print("Done!\r")
    end
end

function time_estimate(dimensions, regimes, neighbourhood_rules, seeds; overwrite=false, path="data/")

    if !isdir(path)
        println("Creating folder...")
        mkdir(path)
    end
    mkPath(distribution_name) = path*distribution_name*"/"

    L = maximum(dimensions)
    t = regimes[1]
    test_seeds = 1:Threads.nthreads()
    
    test_time = @elapsed begin
        
        for neighbourhood_rule in neighbourhood_rules
            distribution_name = "t=$t Uniform" * (neighbourhood_rule=="" ? "" : " " * neighbourhood_rule)
            println("Distribution: $distribution_name, L: $L          ")
            generate_data(mkPath(distribution_name),L, test_seeds, distribution_name, t, overwrite, neighbourhood_rule)
        end
    end
    time_estimate = test_time * length(regimes) * length(seeds)/length(test_seeds)
    formated_time = Dates.canonicalize(Dates.CompoundPeriod(Dates.Second(floor(Int64, time_estimate))))
    println("This will probably take: $formated_time")
end


function itterate_settings(dimensions, regimes, neighbourhood_rules, seeds; overwrite=false, path="data/", estimate_time=true)

    if !isdir(path)
        println("Creating folder...")
        mkdir(path)
    end
    mkPath(distribution_name) = path*distribution_name*"/"

    if estimate_time
        println("Estimating time of run")
        time_estimate(dimensions, regimes, neighbourhood_rules, seeds, overwrite=overwrite, path=path)
    end

    for L in dimensions
        for t in regimes
            for neighbourhood_rule in neighbourhood_rules
                distribution_name = "t=$t Uniform" * (neighbourhood_rule=="" ? "" : " " * neighbourhood_rule)
                println("Distribution: $distribution_name, L: $L          ")
                generate_data(mkPath(distribution_name),L, seeds, distribution_name, t, overwrite, neighbourhood_rule)
            end
        end
    end
end