include("showBundle.jl")

function show_spanning_cluster()
    nr = "SNR"
    path = "data/"
    t = 0.1
    L=1024
    α = 2.0
    seed = 2
    setting = make_settings("Uniform", L, t, nr, α, path)
    file = load_file(setting, average=false)
    b = get_fb(L, nr=nr, without_storage=true)
    b.status = file["spanning_cluster_state/$seed"]

    shift_spanning_cluster!(b)
    resetBundle!(b)
    update_σ!(b)
    p = plot_fb(b, show=false)
    display(p)
end

show_spanning_cluster()