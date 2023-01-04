include("showBundle.jl")
include("../support/inertia.jl")

function show_spanning_cluster()
    nr = "CLS"
    t = 0.9
    L=8
    α = 2.0
    seeds = 10:20
    for seed=seeds, L=[8,512]
        settings = make_settings(L, t, nr, α)
        b = get_bundles_from_settings(settings, seeds=seed, step=-0)
        R = find_radius_of_gyration(b)
        p = plot_fb(b, show=false)
        #minor_axes, major_axes, minor_values, major_values = find_major_and_minor_axes(b)
        #plot_fb_axes(b, minor_axes, major_axes, minor_values, major_values)
        plot_gyration_radi(b, R, nr=1)
        display(p)
    end
end

function show_progression()
    nr = "CLS"
    t = 0.1
    L=256
    α = 1.3
    seed = 1
    settings = make_settings(L, t, nr, α)
    b = get_bundles_from_settings(settings, seeds=seed, step=-0)
    p = plot_fb(b, show=false)


end

show_spanning_cluster()