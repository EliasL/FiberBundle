include("ploting_settings.jl")
include("../burningMan.jl")
include("../support/dataManager.jl")

function forceOnBundle(b::FB, s::FBS)
    for _ in 1:b.N-1
        findAndBreakNextFiber!(b, s)
        update_tension!(b)
        resetBundle!(b)                
    end
    p = plot(1:b.N, s.most_stressed_fiber)
    return p
end

nr = "LLS"
t = 0.0
L=6
α = 2.0
seed = 1
dist="ConstantAverageUniform"
data_path="newData/"

settings = make_settings(l, t, nr, α, dist, data_path)
b,s = get_bundles_from_settings(settings, without_storage=false)
p = forceOnBundle(b,s)
savefig("plots/Graphs/ForceOfSingleBundle")