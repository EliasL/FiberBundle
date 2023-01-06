include("showBundle.jl")
include("../support/inertia.jl")
include("../dataGenerator.jl")

function show_spanning_cluster(L, T, α, nr, seeds)
    for seed=seeds, l=L, t=T
        settings = make_settings(l, t, nr, α)
        b = get_bundles_from_settings(settings, seeds=seed, spanning=true)
        R = find_radius_of_gyration(b)
        p = plot_fb(b, show=false)
        #minor_axes, major_axes, minor_values, major_values = find_major_and_minor_axes(b)
        #plot_fb_axes(b, minor_axes, major_axes, minor_values, major_values)
        plot_gyration_radi(b, R, nr=1)
        display(p)
    end
end

function show_progression(L, t, α, nr, seed)

    steps = 4
    progression = round.((1:steps) ./ (steps+1), digits=2)
    plots = []
    for progress in progression
        settings = make_settings(L, t, nr, α)
        b = get_bundles_from_settings(settings, seeds=seed, progression=progress)
        p = plot_fb(b, show=false)
        title!(p, L"k/N="*"$progress")
        push!(plots, p)
    end
    p = plot(plots...)
    display(p)
    savefig("plots/Visualizations/Progressions/$(L)$(nr)_$(t)_$(α).pdf")  
end

nr = "ELS"
t = [0.22, 0.24]
L=32
α = 2.0
seeds = 1:5
show_spanning_cluster(L, t, α, nr, seeds)
#= s = make_settings(L, t, nr, α)
break_bundle(s, nothing, nothing, seeds, use_threads=false)
show_progression(L, t, α, nr, seeds) =#