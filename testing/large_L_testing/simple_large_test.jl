using Distributed
using ProgressMeter
using Logging
using Profile
using PProf


include("../../support/logingLevels.jl")

@logmsg nodeLog "Preparing workers..."

using JLD2, CodecLz4, Logging
include("../../burningMan.jl")
include("../../support/dataManager.jl")
include("../../support/distributions.jl")


function break_bundle(L, α, distribution::Function, progress_channel, working_channel,
                                  file_name, neighbourhood_rule; seed=0, save_data=true, use_threads=true)
    if use_threads
        put!(working_channel, true) # Indicate a process has started
    end
    N = L*L # Number of fibers
    @assert seed != -1 ""
    Random.seed!(seed)
    x = distribution(N) # Max extension (Distribution)
    neighbours = fillAdjacent(L, NEIGHBOURS)
    neighbourhoods = fillAdjacent(L, NEIGHBOURHOOD)#
    neighbourhood_values = zeros(Int64, N)

    work_values = zeros(Float64, N)
    # These values are reset for each step
    σ  = ones(Float64, N) # Relative tension
    tension = zeros(Float64, N)
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
        i = findNextFiber(tension, σ, x)
        max_σ = tension[i]
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
                tension_storage[storage_index, :] = tension
                if storage_index < length(steps_to_store)
                    storage_index += 1  
                end
            end
        end
        if spanning_cluster != -1 && spanning_cluster_has_not_been_found
            spanning_cluster_state_storage = copy(status)
            spanning_cluster_tension_storage = tension
            spanning_cluster_size_storage = spanning_cluster_size
            spanning_cluster_perimiter_storage = cluster_outline_length[spanning_cluster]
            spanning_cluster_step = step
            spanning_cluster_has_not_been_found = false
        end
    end
end



L=256
α=2.0
t=0.0
dist="Uniform"
nr = "SNR"
path = "data/"
seed=0
distribution = get_uniform_distribution(t)
progress_channel=nothing
working_channel=nothing
file_name = get_file_name(make_settings(dist, L, t, nr, α, path), seed)

@time break_bundle(L, α, distribution, progress_channel, working_channel,
    file_name, nr; seed=0, save_data=true, use_threads=false)