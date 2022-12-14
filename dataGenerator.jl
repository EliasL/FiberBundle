using Distributed
using ProgressMeter
using Logging
using Dates

include("support/logingLevels.jl")

using JLD2, CodecLz4, Logging
include("burningMan.jl")
include("support/dataManager.jl")
include("support/distributions.jl")

function break_bundle(settings, progress_channel, working_channel, seed;
    save_data=true, use_threads=true, stop_after_spanning=false, use_past_progress=false)
    
    if use_threads
        put!(working_channel, true) # Indicate a process has started
    end


    file_name = get_file_name(settings, seed)       

    # Check if there already exists previous data
    if use_past_progress
        @assert false "This feature does not work!"
        b, s = get_bundles_from_settings(settings, seeds=seed, without_storage=false, update_tension=false)
    else
        b, s = get_fb(settings, seed)
    end


    # Break the bundle
    simulation_time = @elapsed for step in b.current_step+1:b.N
        # Simulate step
        findAndBreakNextFiber!(b, s)
        if stop_after_spanning && s.spanning_cluster_has_been_found
            break
        end

        if use_threads
           put!(progress_channel, true) # trigger a progress bar update
        end
    end

    if save_data
        jldopen(file_name, "w") do file
            file["last_step"] = b.current_step
            file["simulation_time"] = simulation_time
            file["spanning_cluster_size"] = s.spanning_cluster_size_storage
            file["spanning_cluster_perimiter"] = s.spanning_cluster_perimiter_storage
            file["spanning_cluster_step"] = s.spanning_cluster_step
            file["most_stressed_fiber"] = s.most_stressed_fiber
            file["nr_clusters"] = s.nr_clusters
            file["break_sequence"] = b.break_sequence#view(b.break_sequence, 1:b.current_step)
            file["largest_cluster"] = s.largest_cluster
            file["largest_perimiter"] = s.largest_perimiter
        end
    end
    if use_threads
        put!(working_channel, false) # trigger a progress bar update
    end
end

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
    flush(stdout)
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
                        flush(stdout)
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

function itterate_settings(dimensions, ??, regimes, neighbourhood_rules, seeds; overwrite=false, path="data/", use_threads=true)
    for L=dimensions, t=regimes, nr=neighbourhood_rules, a=??
        # There is no point in itterating over alphas when using LLS
        settings = make_settings(L, t, nr, a, path)
        @logmsg settingLog "$(now()): Starting $(settings["name"])"
        generate_data(settings, seeds, overwrite; use_threads=use_threads)
    end
end