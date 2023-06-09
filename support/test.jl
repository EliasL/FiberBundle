include("bundleAnalasys.jl")
include("dataManager.jl")
include("../plottingScripts/showBundle.jl")

function do_box_count(l, nr, t, file)
    for seed in 1:1
        b = get_bundle_from_file(file, l, nr=nr, t=t, spanning=true, seed=seed)
        counts = box_counting(b)
        println(counts)
    end
end

L = 512
α = 2.0
t = 0.5
nr = "LLS"
dist = "ConstantAverageUniform"
dataPath = "newData/"

bulk_file = load_file(L, α, t, nr, dist, data_path=dataPath, average=false)
do_box_count(L, nr, t, bulk_file)