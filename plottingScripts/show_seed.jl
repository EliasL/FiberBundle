using LaTeXStrings

include("showBundle.jl")
include("../support/inertia.jl")
include("../dataGenerator.jl")
include("ploting_settings.jl")

function show_spanning_cluster()
    nr = "ELS"
    t = 0.0
    L=128
    α = 2.0
    seeds = 1:5
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

function show_progression()
    NR = ["ELS", "LLS", "CLS"]
    T = 0.0
    L=128
    α = 2.0
    seeds = 0

    for seed=seeds, l=L, t=T, nr=NR
        steps = 4
        progression = round.((1:steps) ./ (steps+1), digits=2)
        plots = []
        for progress in progression
            settings = make_settings(l, t, nr, α)
            b = get_bundles_from_settings(settings, seeds=seed, progression=progress)
            p = plot_fb(b, show=false)
            title!(p, L"k/N="*"$progress")
            push!(plots, p)
        end
        plot(plots..., layout=(1,steps), size=(L*(steps+1), L*1.4))
        savefig("plots/Visualizations/Progressions/$(L)$(nr)_$(t)_$(α).pdf")  
    end
end


function generate_illustrations()
    
    save_path = "plots/Visualizations/differenceIllustrations"
    t = 0.0
    L=128
    α = 2.0
    seed = 1
    nr = "LLS"

    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed", save_path)
    nr="CLS"
    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed", save_path)
    L=1024
    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed", save_path)
    nr="LLS"
    save_picture(L, nr, t, α, seed, "$(nr)$(L)s$seed", save_path)
end

generate_illustrations()

#show_spanning_cluster()
#show_progression()