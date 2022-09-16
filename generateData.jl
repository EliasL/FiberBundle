
using JLD2
include("burningMan.jl")

function get_name(L, distribution, seed)
    return path*distribution*", size $L x $L, seed: $seed.jld"
end

function save_step(file, step, x, σ, status, cluster_size, cluster_outline, cluster_outline_length)
    group = JLD2.Group(file, "$step")
    group["x"] = x
    group["σ"] = σ
    group["status"] = status
    group["cluster_size"] = cluster_size
    group["cluster_outline"] = cluster_outline
    group["cluster_outline_length"] = cluster_outline_length
end

function main(L; distribution = "Uniform", seed::Int64=0, show_progress=true)

    # First check if we have already run this simulation
    name = get_name(L, distribution, seed[])
    if isfile(name) && false
        println("\nFound $name")
        return
    end


    # We will continously be writing to our file in order to save RAM
    jldopen(name, "w"; compress=false) do file

        N = L*L # Number of fibers
        Random.seed!(seed)
        x = rand(Float64, N) # Max extension (Distribution)
        σ  = ones(Float64, N) # Relative tension
        adjacent = fillAdjacent(L)
        status = fill(-1,N)
        cluster_size = spzeros(Int64, N)
        cluster_outline = spzeros(Int64, N)
        cluster_outline_length = spzeros(Int64, N)
        unexplored = zeros(Int64, N)

        # Do what you want
        p = Progress(N; enabled=show_progress)
        for step in 1:N

            i = findNextFiber(σ, x)
            resetClusters(status, σ)
            break_fiber(i, status, σ)
            update_σ(status,σ,adjacent, cluster_size, cluster_outline, cluster_outline_length, unexplored)

            # Save data
            # We don't want to wait for the pc to write the data, so we copy the data and we make a 
            # subtask to handle the writing to file. 
            # Here there is a balance of cores vs ram. If you have loads of ram, just run multiple instances
            # if you have loads of cores, you can use the spare cores to multithread data saving
            Threads.@spawn save_step(file, step, map(copy, [x, σ, status, cluster_size, cluster_outline, cluster_outline_length])...)

            next!(p)

        end #step loop
    end #file
end

function worker(seeds::Vector{Int64})
    for s in seeds
        main(L; seed=s, show_progress=false)
    end
end 

const L = 64
const threads = 3
const seeds = 1000
const progress = Progress(seeds)
const path = "data/"

chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]

println(Threads.nthreads())
for seeds in chunk(1:seeds, threads)
    Threads.@spawn main(L; seed=seeds, show_progress=true)
end
