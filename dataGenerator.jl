using Distributed
using ProgressMeter
using Logging

include("support/logingLevels.jl")

@logmsg nodeLog "Preparing workers..."

using JLD2, CodecLz4, Logging
include("burningMan.jl")
include("support/dataManager.jl")
include("support/distributions.jl")

function break_bundle(settings, progress_channel, working_channel, seed;
    save_data=true, use_threads=true, stop_after_spanning=true)
    
    if use_threads
        put!(working_channel, true) # Indicate a process has started
    end


    file_name = get_file_name(settings, seed)       
    @assert seed != -1 ""
    Random.seed!(seed)
    
    b::FB, s::FBS = get_fb(settings)
    # Break the bundle
    @time simulation_time = @elapsed for step in 1:b.N
        # Simulate step
        findNextFiber!(b)
        resetBundle!(b)
        break_fiber!(b)
        update_σ!(b)
        
        update_storage!(b, s, seed)
        if stop_after_spanning && s.spanning_cluster_has_been_found
            break
        end

        if use_threads
           put!(progress_channel, true) # trigger a progress bar update
        end
    end

    if save_data
        jldopen(file_name, "w") do file
            if seed <= 10
                file["sample_states"] = s.status_storage
                file["tension"] = s.tension_storage
                file["spanning_cluster_state"] = s.spanning_cluster_state_storage
                file["spanning_cluster_tension"] = s.spanning_cluster_tension_storage
            end
            file["last_step"] = b.current_step
            file["simulation_time"] = simulation_time
            file["spanning_cluster_size"] = s.spanning_cluster_size_storage
            file["spanning_cluster_perimiter"] = s.spanning_cluster_perimiter_storage
            file["spanning_cluster_step"] = s.spanning_cluster_step
            file["sample_states_steps"] = s.steps_to_store
            file["most_stressed_fiber"] = s.most_stressed_fiber
            file["nr_clusters"] = s.nr_clusters
            file["break_sequence"] = b.break_sequence
            file["largest_cluster"] = s.largest_cluster
            file["largest_perimiter"] = s.largest_perimiter
        end
    end
    if use_threads
        put!(working_channel, false) # trigger a progress bar update
    end
end

# Done preparing workers!
@logmsg nodeLog "Done!"


function run_workers(settings, seeds; save_data=true, use_threads=true)
    
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

    if use_threads
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
        
            # the second task does the computation
            @async begin
                @distributed (+) for i in seeds
                break_bundle(settings, progress, working, i; save_data=save_data)
                i^2 #I have no idea what this does
                end
            end
end             
    else
        @showprogress for i in seeds
            break_bundle(settings, progress, working, i; save_data=save_data, use_threads=false)
        end
    end
end

function generate_data(settings, requested_seeds, overwrite; save_data=true, use_threads=true)
    # If we don't want to save the data, we just do this
    if !save_data
        run_workers(settings, requested_seeds, save_data=save_data, use_threads=use_threads)
        return
    end

    missing_seeds = prepare_run(settings, requested_seeds, overwrite)

    if length(missing_seeds) > 0
        @logmsg settingLog "Running workers..."
        run_workers(settings, missing_seeds, save_data=save_data, use_threads=use_threads)
        @logmsg settingLog "Done!"


        @logmsg settingLog "Calculating averages... "
        clean_after_run(settings, requested_seeds)
        @logmsg settingLog "Done!"
    end
end

function itterate_settings(dimensions, α, regimes, neighbourhood_rules, seeds; overwrite=false, path="data/", use_threads=true)
    for L=dimensions, t=regimes, nr=neighbourhood_rules, a=α
        # There is no point in itterating over alphas when using UNR
        if nr=="UNR"
            a = 0.0
        end
        settings = make_settings("Uniform", L, t, nr, a, path)
        @logmsg settingLog "Starting $(settings["name"])"
        generate_data(settings, seeds, overwrite; use_threads=use_threads)
    end
end