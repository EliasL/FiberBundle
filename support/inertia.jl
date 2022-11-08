include("dataManager.jl")
include("../burningMan.jl")
include("../plottingScripts/showBundle.jl")



function find_center_of_mass(b::FB)
    update_σ!(b)
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
    #println(file["spanning_cluster_state/$seed"])
    b.status = file["spanning_cluster_state/$seed"]
    update_σ!(b)
    show_fb(b)
end

test()