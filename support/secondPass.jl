#Some properties might be interesting, but we chose not to store it

include("../plottingScripts/ploting_settings.jl")
include("../burningMan.jl")
include("../support/dataManager.jl")

function forceOnBundle(b::FB, s::FBS)
    x = zeros(Float64, b.N)
    σ = zeros(Float64, b.N)
    for i in 1:b.N
        update_tension!(b)
        find_next_fiber!(b)
        update_storage!(b, s)
        f = b.break_sequence[i]
        X = b.x[f]#/b.σ[f] 
        x[i] = i
        #x[2*i+1] = x[2*i]
        σ[i] = (b.N-i) *b.max_σ
        #σ[2*i+1] = (b.N-i-1) * b.max_σ
        break_fiber!(b)
        
        resetBundle!(b)
    end
    p = plot(x,σ, legend=:topleft, title=b.nr)
    return p
end


function load_data(bulk_file, key, l, nr, t, dist, xKey, xKeyf,data_path)
    path = "$data_path/$xKey/"
    name = "$(dist)_$(l)_$(nr)_$(t)_$(key)_$(string(xKeyf))"
    avg_x = 0
    if !isdir(path)
        mkpath(path)
    end
    if isfile("$(path)$(name).jld2")
        f = load("$(path)$(name).jld2")
        value = f[name]
        avg_x = f[name*"x"]
    else
        value = 0
        seeds_used = bulk_file["seeds_used"]
        nr_seeds = length(seeds_used)
        for seed in bulk_file["seeds_used"]
            x = xKeyf(bulk_file["$xKey/$seed"])
            value += bulk_file["$key/$seed"][x]
            avg_x += x
        end
        value /= nr_seeds
        avg_x /= nr_seeds
        # Save data
        jldopen("$(path)$(name).jld2", "w") do file
            file[name] = value
            file[name*"x"] = avg_x
        end
    end

    return value, avg_x
end

nr = "ELS"
t = 0.5
L=6
α = 2.0
seed = 1
dist="ConstantAverageUniform"
data_path="newData/"