using Distributed
using ProgressMeter

rmprocs(workers())
threads = Threads.nthreads()

print("Starting $threads workers... ")
addprocs(threads)
println("Done!")
@everywhere begin
    using CodecBzip2
    using JLD2
    using Pkg
end


@everywhere include("burningMan.jl")
@everywhere include("dataManager.jl")
@everywhere include("distributions.jl")

@everywhere function break_bundle(L, distribution::Function, progress_channel, working_channel, file_name, using_neighbourhood_rules; seed=0)
    put!(working_channel, true) # Indicate a process has started

    N = L*L # Number of fibers
    @assert seed != -1 ""
    Random.seed!(seed)
    x = distribution(N) # Max extension (Distribution)
    neighbours = fillAdjacent(L, NEIGHBOURS)
    neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)

    # These values are reset for each step
    σ  = ones(Float64, N) # Relative tension
    status = fill(-1,N)
    cluster_size = zeros(Int64, N)
    cluster_outline_length = zeros(Int64, N)
    # These values are reset for each cluster
    cluster_outline = zeros(Int64, N)
    unexplored = zeros(Int64, N)

    # These arrays store one value for each step
    steps = N # just for readability
    most_stressed_fiber = zeros(Int64, steps)
    nr_clusters = zeros(Int64, steps)
    largest_cluster = zeros(Int64, steps)
    largest_perimiter = zeros(Int64, steps)

    # We want to store some samples of the processed
    # I'm thinking at 10%, 20%, ... 90% done would work
    # ie, 9 images
    division = 10
    status_storage = zeros(Int64, division-1, N)
    # If N=100 Steps to store is now [90, 80, ... , 10]
    steps_to_store = [round(Int64,N/division * i) for i in 1:division-1]
    storage_index = 1

    # Do what you want
    for step in 1:N

        # Simulate step
        i = findNextFiber(σ, x)
        resetClusters(status, σ)
        break_fiber(i, status, σ)
        _nr_clusters = update_σ(status, σ, neighbours, neighbourhoods, cluster_size, cluster_outline, cluster_outline_length, unexplored; use_neighbourhood_rules=using_neighbourhood_rules)

        # Save important data from step
        most_stressed_fiber[step] = σ[i]/x[i]
        nr_clusters[step] = _nr_clusters
        largest_cluster[step] = maximum(cluster_size)
        largest_perimiter[step] = maximum(cluster_outline_length)

        # Save step for visualization
        if step == steps_to_store[storage_index]
            status_storage[storage_index, :] = status
            if storage_index < length(steps_to_store)
                storage_index += 1  
            end
        end

        put!(progress_channel, true) # trigger a progress bar update
    end 
    
    jldopen(file_name, "w") do file
        file["sample_states"] = status_storage
        file["sample_states_steps"] = steps_to_store./N
        file["most_stressed_fiber"] = most_stressed_fiber./N
        file["nr_clusters"] = nr_clusters./N
        file["largest_cluster"] = largest_cluster./N
        file["largest_perimiter"] = largest_perimiter./N
    end
    put!(working_channel, false) # trigger a progress bar update

end

function run_workers(L, distribution_name, distribution_function, seeds, path, using_neighbourhood_rules)
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

        # the second task does the computation
        @async begin
            @distributed (+) for i in seeds
                name = get_name(L, distribution_name, path, i)
                break_bundle(L, distribution_function, progress, working, name, using_neighbourhood_rules; seed=i)
                i^2
            end
        end
    end
end

function generate_data(path, L, requested_seeds, distribution_name, overwrite, using_neighbourhood_rules)

    println("Processing L = $L ...")

    if !isdir(path)
        println("Creating folder... ")
        mkdir(path)
    end

    missing_seeds = prepare_run(L, distribution_name, path, requested_seeds, overwrite)

    # get distribtion function
    distribution_function(n) = zeros(n) # Default

    if distribution_name == "Uniform with Neighbourhood rules"
        distribution_function = uniform
    end

    if length(missing_seeds) > 0
        println("Running workers...")
        run_workers(L, distribution_name, distribution_function, missing_seeds, path, using_neighbourhood_rules)
        println("Done!")


        println("Calculating averages ...")
        clean_after_run(L, name_function(), name_function, requested_seeds)
        println("Done!")
    end
end


seeds = 1:200
distribution_name = "Uniform with Neighbourhood rules"
using_neighbourhood_rules = true
overwrite = false
global_path = "data/"
if !isdir(global_path)
    println("Creating folder...")
    mkdir(global_path)
end
mkPath(L) = global_path*distribution_name*"/"

for L in [256,32, 64, 128] 
    generate_data(mkPath(L),L, seeds, distribution_name, overwrite, using_neighbourhood_rules)
end



println("Removing workers")
rmprocs(workers())