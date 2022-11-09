include("dataManager.jl")
include("../burningMan.jl")
include("../plottingScripts/showBundle.jl")


function find_center_of_mass(b::FB, only_spanning=true)
    # NB THIS ONLY WORKS FOR CLUSTERS THAT ARE NOT ON THE PERIODIC BORDER
    
    avg_i = zeros(Float64, b.N)
    avg_j = zeros(Float64, b.N)
    
    
    resetBundle!(b)
    # Redefinition of update_σ
    # For every fiber in the plane
    for fiber in eachindex(b.status)
        #@logmsg σUpdateLog "Breaking fiber $i"
        # If it is broken and unexplored
        if b.status[fiber] == 0#BROKEN
            # We have found a new cluster!
            # increase number of clusters
            b.c += 1
            # assign fiber to this new cluster
            b.status[fiber] = b.c
            # explore the new cluster
            #@logmsg σUpdateLog "Exploring cluster"
            found_spanning_cluster = explore_cluster_at!(fiber, b)
            if found_spanning_cluster && b.spanning_cluster_id == -1
                b.spanning_cluster_id = b.c
            end

            # We should now have updated cluster_outline,
            # and with that we can update sigma for one cluster
            update_cluster_outline_stress!(b)


            # index to coordinates
            y = mod1(fiber,b.L)
            x = floor(Int64, fiber/b.L)

            avg_i[b.c] +=  y
            avg_j[b.c] +=  x
            
            for i in 1:b.cluster_size[b.c]
                
                avg_i[b.c] += b.rel_pos_y[i] + y
                avg_j[b.c] += b.rel_pos_x[i] + x
            end
            avg_i[b.c] /= b.cluster_size[b.c]
            avg_j[b.c] /= b.cluster_size[b.c]
        end

    end
    #println(b.spanning_cluster_id)
    #println(avg_j[b.spanning_cluster_id], ", ", avg_i[b.spanning_cluster_id])
    if only_spanning
        return avg_j[b.spanning_cluster_id], avg_i[b.spanning_cluster_id]
    else
        return avg_j[1:b.c], avg_i[1:b.c]
    end
end

function find_center_of_mass1(b::FB, only_spanning=true)
    # NB THIS ONLY WORKS FOR CLUSTERS THAT ARE NOT ON THE PERIODIC BORDER
    resetBundle!(b)
    update_σ!(b)
    nr_clusters = b.c
    avg_i = zeros(Float64, nr_clusters)
    avg_j = zeros(Float64, nr_clusters)
    test_x = []
    test_y = []
    weights = 1 ./ b.cluster_size
    for fiber in eachindex(b.status)
        
        # index to coordinates
        y = mod1(fiber, b.L)
        x = floor(Int64, fiber/b.L)
        rel_x = b.rel_pos_x[fiber]
        rel_y = b.rel_pos_y[fiber]
        # Check if the fiber is part of a cluster or not
        cluster = b.status[fiber]
        if cluster > 0 && (cluster == b.spanning_cluster_id || !only_spanning)
            push!(test_x,x + rel_x)
            push!(test_y,y + rel_y)
            if rel_x == 0 && rel_y == 0
                avg_i[cluster] += mod1(fiber, b.L)
                avg_j[cluster] += floor(Int64, fiber/b.L)
            else
                avg_i[cluster] += (rel_x) * weights[cluster]
                avg_j[cluster] += (rel_y) * weights[cluster]
            end
        end
    end

    return test_x, test_y#avg_j, avg_i
end

function test()
    nr = "SNR"
    path = "data/"
    t = 0.1
    L=256
    α = 2.0
    seed = 1
    setting = make_settings("Uniform", L, t, nr, α, path)
    file = load_file(setting, average=false)
    b = get_fb(L, nr=nr, without_storage=true)
    b.status = file["spanning_cluster_state/$seed"]
    shift_spanning_cluster!(b)
    show_fb(b)
    cm = find_center_of_mass(b, true)
    plot!(cm, seriestype = :scatter)
end

test()