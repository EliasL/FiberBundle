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
    L = settings["L"]
    α = settings["a"]
    t = settings["t"]
    nr = settings["nr"]
    dist = settings["dist"]
    file_name = get_file_name(settings, seed)

    @assert seed != -1 ""
    Random.seed!(seed)

    b, s = get_fb(L,α=α, t=t, nr=nr, dist=dist)
    # Break the bundle
    simulation_time = @elapsed for step in 1:b.N
        # Simulate step
            if false
                @time findNextFiber!(b)
                @time resetBundle!(b)
                @time break_fiber!(b)
                @time update_σ!(b)
                @time update_storage!(b, s)
            else 
                findNextFiber!(b)
                resetBundle!(b)
                break_fiber!(b)
                update_σ!(b)
                update_storage!(b, s)
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



L=32
α=2.0
t=0.0
nr = "CLS"
seed=0
settings = make_settings(L, t, nr, α)

@time break_bundle(settings, seed, save_data=true)
println("ok")