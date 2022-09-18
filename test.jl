using ProgressMeter
using Distributed

@everywhere using JLD2, CodecBzip2


@everywhere include("burningMan.jl")


@everywhere function get_name(L, distribution, seed::Int64)
    #return path*distribution*", size $L x $L, seed: $seed.jld"
    return path*distribution*", size $L x $L, seed $seed.jld"
end

@everywhere function save_step(file, step, x, σ, status)
    file["x/$step"] = copy(x)
    file["σ/$step"] = copy(σ)
    file["status/$step"] = copy(status)
end

@everywhere function main(L; distribution = "Uniform", seed=0, show_progress=true)
    # First check if we have already run this simulation
    name = get_name(L, distribution, seed)
    if !isdir(path)
        println("Creating folder...")
        mkdir(path)
    end
    if isfile(name) && false
        println("\nFound $name")
        return
    end
    # We will continously be writing to our file in order to save RAM
    # use false or Bzip2Compressor(). This saves 35% space but increases runtime by 700%
    # so don't use it
    
    jldopen(name, "w"; compress=false) do file
        N = L*L # Number of fibers
        Random.seed!(seed)
        x = rand(Float64, N) # Max extension (Distribution)
        σ  = ones(Float64, N) # Relative tension
        adjacent = fillAdjacent(L)
        status = fill(-1,N)
        cluster_size = spzeros(Int64, N, N)
        cluster_outline = spzeros(Int64, N, N)
        cluster_outline_length = spzeros(Int64, N, N)
        unexplored = zeros(Int64, N)

        # Do what you want
        p = Progress(N; enabled=show_progress)
        for step in 1:N

            i = findNextFiber(σ, x)
            resetClusters(status, σ)
            break_fiber(i, status, σ)
            update_σ(status,σ,adjacent, cluster_size[step, :], cluster_outline[step, :], cluster_outline_length[step, :], unexplored)

            # Save data
            # We don't want to wait for the pc to write the data, so we copy the data and we make a 
            # subtask to handle the writing to file. 
            # Here there is a balance of cores vs ram. If you have loads of ram, just run multiple instances
            # if you have loads of cores, you can use the spare cores to multithread data saving
            save_step(file, step, x, σ, status)
            next!(p)
        end #step loop
        # Save the rest of the data
        file["cluster_size"] = cluster_size
        file["cluster_outline"] = cluster_outline
        file["cluster_outline_length"] = cluster_outline_length
    end #file
end
@everywhere path = "data/"
n_steps = 20
p = Progress(n_steps)
channel = RemoteChannel(()->Channel{Bool}(), 1)
threads = Threads.nthreads()#Int(ceil(Threads.nthreads()/2))
addprocs(threads-1)
println("Threads: $threads")

@sync begin # start two tasks which will be synced in the very end
    # the first task updates the progress bar
    @async while take!(channel)
        next!(p)
    end

    # the second task does the computation
    @async begin
        @distributed (+) for i in 1:n_steps
            main(15, seed=i, show_progress=false)
            put!(channel, true) # trigger a progress bar update
            i^2
        end
        put!(channel, false) # this tells the printing task to finish
    end
end