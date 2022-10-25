using Distributed
using ProgressMeter
using Logging

include("support/logingLevels.jl")

@logmsg nodeLog "Preparing workers..."

@everywhere using JLD2, CodecLz4, Logging
@everywhere include("burningMan.jl")
@everywhere include("support/dataManager.jl")
@everywhere include("support/distributions.jl")

@everywhere function break_bundle(L, α, distribution::Function, progress_channel, working_channel,
                                  file_name, neighbourhood_rule; seed=0, save_data=true)
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
    # Relative possition of every fiber with respect to it's cluster
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
    spanning_cluster_state_storage = zeros(Int64, N)
    spanning_cluster_tension_storage = zeros(Int64, N)
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

    @logmsg threadLog "Starting seed $seed..."
    for step in 1:N
        # Simulate step
        i = findNextFiber(σ, x)
        max_σ = σ[i]/x[i]
        resetClusters(status, σ)
        @logmsg singleBreakLog "Breaking fiber nr $step"
        break_fiber(i, status, σ)
        @logmsg singleBreakLog "Updating cluster"
        _nr_clusters, spanning_cluster, spanning_cluster_size = update_σ(status, σ, neighbours, neighbourhoods,
                                                                         neighbourhood_values, cluster_size,
                                                                         cluster_dimensions, rel_pos_x, rel_pos_y,
                                                                         cluster_outline, cluster_outline_length,
                                                                         unexplored; neighbourhood_rule=neighbourhood_rule, α=α)
        
        @logmsg singleBreakLog "Save data"
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
            spanning_cluster_state_storage = copy(status)
            spanning_cluster_tension_storage = σ ./ x
            spanning_cluster_size_storage = spanning_cluster_size
            spanning_cluster_perimiter_storage = cluster_outline_length[spanning_cluster]
            spanning_cluster_step = step
            spanning_cluster_has_not_been_found = false
        end
        @debug "adding progress"
        put!(progress_channel, true) # trigger a progress bar update
        @debug "step is done"
    end
    @logmsg threadLog "Saving data..."
    if save_data
        jldopen(file_name, "w") do file
            if seed <= 10
                file["sample_states"] = status_storage
                file["tension"] = tension_storage
                file["spanning_cluster_state"] = spanning_cluster_state_storage
                file["spanning_cluster_tension"] = spanning_cluster_tension_storage
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
    end
    put!(working_channel, false) # trigger a progress bar update

end

# Done preparing workers!
@logmsg nodeLog "Done!"


function run_workers(settings, distribution_function, seeds; save_data=true)
    
    # Check if the current logger level is set to log on the thread level
    L = settings["L"]
    show_progress = Logging.min_enabled_level(global_logger()) == threadLog
    p = Progress(length(seeds)*L^2, enabled=show_progress)
    
    p = 0
    print_step = 0
    finish = length(seeds)*L^2

    progress = RemoteChannel(()->Channel{Bool}(), 1)
    working = RemoteChannel(()->Channel{Bool}(), 1)
    active_workers = Threads.Atomic{Int}(0)
    completed_runs = Threads.Atomic{Int}(0)

    @logmsg settingLog "Starting work on $(settings["name"])"
    @sync begin # start two tasks which will be synced in the very end
        # the first task updates the progress bar
        @async begin
            while take!(progress)
                #ProgressMeter.next!(p; showvalues = [("Active workers", active_workers[]), ("Completed tasks", completed_runs[])])
                p += 1
                percent_done = p/finish*100
                if percent_done > print_step
                    print_step += 1
                    @logmsg settingLog "$(print_step-1)%, Active workers: $(active_workers[]), Completed tasks: $(completed_runs[])."
                end
            end
            @logmsg settingLog "100%, Active workers: $(active_workers[]), Completed tasks: $(completed_runs[])."
            #ProgressMeter.next!(p; showvalues = [("Active workers", active_workers[]), ("Completed tasks", completed_runs[])])
            #ProgressMeter.finish!(p; showvalues = [("Active workers", active_workers[]), ("Completed tasks", completed_runs[])])
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
                name = get_file_name(settings, i)
                break_bundle(settings["L"], settings["a"], distribution_function, progress, working, name, settings["nr"]; seed=i, save_data=save_data)
                i^2 #I have no idea what this does
            end
        end
    end
end

function generate_data(settings, requested_seeds, overwrite; save_data=true)
    # get distribtion function
    distribution_function(n) = nothing

    if  occursin("Uniform", settings["name"])
        distribution_function = get_uniform_distribution(settings["t"])
    else
        error("No distribution found! Got: $distribution_name")
    end

    # If we don't want to save the data, we just do this
    if !save_data
        run_workers(settings, distribution_function, requested_seeds, save_data=save_data)
        return
    end

    missing_seeds = prepare_run(settings, requested_seeds, overwrite)

    if length(missing_seeds) > 0
        @logmsg settingLog "Running workers..."
        run_workers(settings, distribution_function, missing_seeds, save_data=save_data)
        @logmsg settingLog "Done!"


        @logmsg settingLog "Calculating averages... "
        clean_after_run(settings, requested_seeds)
        @logmsg settingLog "Done!"
    end
end

function itterate_settings(dimensions, α, regimes, neighbourhood_rules, seeds; overwrite=false, path="data/")
    for L=dimensions, t=regimes, nr=neighbourhood_rules, a=α

        settings = make_settings("Uniform", L, t, nr, a, path)
        @logmsg settingLog "Starting $(settings["name"])"
        generate_data(settings, seeds, overwrite)
    end
end