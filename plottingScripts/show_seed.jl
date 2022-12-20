include("showBundle.jl")
include("../support/inertia.jl")

function show_spanning_cluster()
    nr = "CLS"
    path = "data/"
    t = 0.1
    L=256
    α = 1.3
    seed = 1
    settings = make_settings("Uniform", L, t, nr, α, path)
    b = get_bundles_from_settings(settings, seeds=seed, step=-0)
    p = plot_fb(b, show=false)
    minor_axes, major_axes, minor_values, major_values = find_major_and_minor_axes(b)
    plot_fb_axes(b, minor_axes, major_axes, minor_values, major_values)
    R = find_radius_of_gyration(b)
    plot_gyration_radi(b, R, nr=2)
    display(p)
end

show_spanning_cluster()