
using JLD2
using Distributed

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

function worker(seed::RemoteChannel, progress::RemoteChannel)
    while true
        main(L; seed=take!(seed), show_progress=false)
        put!(progress, true)
    end

end 

L = 64
threads = 3
seeds = 1000
p = Progress(seeds)
progress_channel = RemoteChannel(()->Channel{Bool}(), 1)
seed_channel = RemoteChannel(()->Channel{Int64}(), 1)
path = "data/"


println(Threads.nthreads())
    
for _ in 1:threads
    Threads.@spawn worker(seed_channel, progress_channel)
end

@sync begin

    println("Start")
    @async for s in 1:seeds
        println("seed")
        put!(seed_channel, s)
    end

    @async while take!(progress_channel)
        println("ok")
        next!(p)
    end
end
