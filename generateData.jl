using Distributed
using ProgressMeter
threads = Threads.nthreads()#Int(ceil(Threads.nthreads()/2))

print("Starting $threads workers... ")
addprocs(threads)
println("Done!")
@everywhere begin
    using CodecBzip2
    using JLD2
    using Pkg
end


@everywhere include("burningMan.jl")


@everywhere function get_name(L, distribution, path, seed=-1)
    if seed == -1
        return path*distribution*"$L.jld2"
    else
        return path*distribution*"$L,$seed.jld2"
    end
end

@everywhere function break_bundle(L, progress_channel, working_channel, file_name; seed=0)
    put!(working_channel, true) # trigger a progress bar update

    N = L*L # Number of fibers
    Random.seed!(seed)
    x = rand(Float64, N) # Max extension (Distribution)
    σ  = ones(Float64, N) # Relative tension
    adjacent = fillAdjacent(L)
    status = fill(-1,N)
    cluster_size = zeros(Int64, N)
    cluster_outline_length = zeros(Int64, N)
    cluster_outline = zeros(Int64, N)
    unexplored = zeros(Int64, N)

    # These arrays store one value for each step
    steps = N
    most_stressed_fiber = zeros(Int64, steps)
    nr_clusters = zeros(Int64, steps)
    largest_cluster = zeros(Int64, steps)
    largest_perimiter = zeros(Int64, steps)

    # Do what you want
    for step in 1:N

        # Simulate step
        i = findNextFiber(σ, x)
        resetClusters(status, σ)
        break_fiber(i, status, σ)
        _nr_clusters = update_σ(status,σ,adjacent, cluster_size, cluster_outline, cluster_outline_length, unexplored)

        # Save important data from step
        most_stressed_fiber[step] = σ[i]/x[i]
        nr_clusters[step] = _nr_clusters
        largest_cluster[step] = maximum(cluster_size)
        largest_perimiter[step] = maximum(cluster_outline_length)


        put!(progress_channel, true) # trigger a progress bar update
    end 
    
    jldopen(file_name, "w") do file
        file["most_stressed_fiber"] = most_stressed_fiber./N
        file["nr_clusters"] = nr_clusters./N
        file["largest_cluster"] = largest_cluster./N
        file["largest_perimiter"] = largest_perimiter./N
    end
    put!(working_channel, false) # trigger a progress bar update

end

function run_workers(L, distribution, seeds, path)
    p = Progress(length(seeds)*L^2)
    progress = RemoteChannel(()->Channel{Bool}(), 1)
    working = RemoteChannel(()->Channel{Bool}(), 1)
    active_workers = Threads.Atomic{Int}(0)
    completed_runs = Threads.Atomic{Int}(nr_seeds-length(seeds))

    @sync begin # start two tasks which will be synced in the very end
        # the first task updates the progress bar
        @async while take!(progress)
            ProgressMeter.next!(p; showvalues = [("Active workers", active_workers[]), ("Completed tasks", completed_runs[])])
        end
        @async while true
            w = take!(working)
            if w
                active_workers[] += 1
            else 
                active_workers[] -=1 
                completed_runs[] += 1
                if nr_seeds == completed_runs[]
                    put!(progress, false) # this tells the printing task to finish
                    break
                end
            end
        end

        # the second task does the computation
        @async begin
            @distributed (+) for i in seeds
                name = get_name(L, distribution, path, i)
                break_bundle(L, progress, working, name, seed=i)
                i^2
            end
        end
    end
end

L = 64
threads = 3
seeds = 1000
p = Progress(seeds)
progress_channel = RemoteChannel(()->Channel{Bool}(), 1)
seed_channel = RemoteChannel(()->Channel{Int64}(), 1)
path = "data/"


function generate_data(path, L, nr_seeds, distribution)

    println("Processing L = $L ...")

    if !isdir(path)
        println("Creating folder... ")
        mkdir(path)
    end

    test_seeds = 1:nr_seeds
    seeds = collect(1:nr_seeds)
    found_files = false
    work_to_do = true
    print("Checking existing seeds... ")

    name = get_name(L, distribution, path)
    if isfile(name) && !overwrite
        jldopen(get_name(L, distribution, path), "r") do file
            for seed in test_seeds
                if haskey(file, "nr_clusters/$seed")
                    if !found_files
                        println("\nFound existing seeds!")
                    end
                    # Write out the file
                    jldopen(get_name(L, distribution, path, seed), "w") do s_file
                        for key in ["most_stressed_fiber", "nr_clusters", "largest_cluster", "largest_perimiter"]
                            s_file[key] = file["$key/$seed"]
                        end
                    end
                    found_files=true
                    #print("$seed, ")
                    deleteat!(seeds, findfirst(x->x==seed,seeds))
                end
            end
        end
    end
    if found_files
        println("\u1b[1D\u1b[1D.")
        if length(seeds) == 0
            println("All seeds have already been processed")
            work_to_do = false
        else
            #println("Processing seeds: " * join(seeds, ", "))
        end
    else
        println("Done!")
    end

    if work_to_do
        println("Running workers...")
        run_workers(L, distribution, seeds, path)
        println("Done!")
    end

    println("Calculate averages ...")
    keys = ["most_stressed_fiber", "nr_clusters", "largest_cluster", "largest_perimiter"]
    averages = Dict()
    for key in keys
        averages[key] = zeros(L*L)
    end
    jldopen(get_name(L, distribution, path), "w") do file
        for seed in 1:nr_seeds
            seed_file_name = get_name(L, distribution, path, seed)
            jldopen(seed_file_name, "r") do s_file
                for key in keys
                    averages[key] += s_file[key] ./ nr_seeds
                    file["$key/$seed"] = s_file[key]
                end
            end
            rm(seed_file_name)
        end
        for key in keys
            file["average_$key"] = averages[key]
        end
        file["seeds_used"] = nr_seeds
    end
    println("Done!")

end


nr_seeds = 500
distribution = "Uniform"
overwrite = false
global_path = "data/"
if !isdir(global_path)
    println("Creating folder...")
    mkdir(global_path)
end
mkPath(L) = global_path*distribution*"/"

for L in [256, 32, 64, 128] 
    generate_data(mkPath(L),L, nr_seeds, distribution)
end



println("Removing workers")
rmprocs(workers())