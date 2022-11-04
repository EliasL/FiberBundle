using Distributed
using ProgressMeter
using Logging
using Profile

include("../../support/logingLevels.jl")

@logmsg nodeLog "Preparing workers..."

using JLD2, CodecLz4, Logging
include("../../burningMan.jl")
include("../../support/dataManager.jl")
include("../../support/distributions.jl")

function break_bundle(settings, seed; save_data=true)
    println("ok1")
    L = settings["L"]
    α = settings["a"]
    t = settings["t"]
    nr = settings["nr"]
    dist = settings["dist"]
    println("ok2")
    file_name = get_file_name(settings, seed)

    N = L*L # Number of fibers
    @assert seed != -1 ""
    Random.seed!(seed)
    println("ok3")

    b, s = get_fb(L,α=α, t=t, nr=nr, dist=dist)
    println("ok4")
    # Break the bundle
    simulation_time = @elapsed @showprogress for step in 1:N
        # Simulate step
        i = findNextFiber!(b)
        resetBundle!(b)
        break_fiber!(i, b)
        update_σ!(b)
        print(step)
        # Save important data from step
        s.most_stressed_fiber[step] = 1/b.max_σ
        s.nr_clusters[step] = b.c # The last cluster id is also the number of clusters
        s.largest_cluster[step] = maximum(b.cluster_size)
        s.largest_perimiter[step] = maximum(b.cluster_outline_length)
        # Save step for visualization
        if seed <= 10 #Only save samples from 10 first seeds
            if step == s.steps_to_store[s.storage_index] &&
            s.storage_index < length(s.steps_to_store)
                s.status_storage[s.storage_index, :] = b.status
                s.tension_storage[s.storage_index, :] = b.tension
                s.storage_index += 1  
            end
        end
        if b.spanning_cluster_id != -1 && s.spanning_cluster_has_not_been_found
            s.spanning_cluster_state_storage = copy(b.status)
            s.spanning_cluster_tension_storage = b.tension
            s.spanning_cluster_size_storage = b.cluster_size[b.spanning_cluster_id]
            s.spanning_cluster_perimiter_storage = b.cluster_outline_length[b.spanning_cluster_id]
            s.spanning_cluster_step = step
            s.spanning_cluster_has_not_been_found = false
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
            file["simulation_time"] = simulation_time
            file["spanning_cluster_size"] = s.spanning_cluster_size_storage
            file["spanning_cluster_perimiter"] = s.spanning_cluster_perimiter_storage
            file["spanning_cluster_step"] = s.spanning_cluster_step
            file["sample_states_steps"] = s.steps_to_store
            file["most_stressed_fiber"] = s.most_stressed_fiber
            file["nr_clusters"] = s.nr_clusters
            file["largest_cluster"] = s.largest_cluster
            file["largest_perimiter"] = s.largest_perimiter
        end
    end
end



L=64
α=2.0
t=0.0
dist="Uniform"
nr = "SNR"
path = "data/"
seed=0
settings = make_settings(dist, L, t, nr, α, path)

println("start")
break_bundle(settings, seed, save_data=true)