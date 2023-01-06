using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")
include("../support/dataManager.jl")


function basicPropertiesPlot(L, ts, nr; use_y_lable=true)
    
    N = L.*L
    files_and_t = []
    for t in ts
        push!(files_and_t, load_file(L, α, t, nr))
        # Check that the data we use is from completely borken bundles
        min_steps = get_min_steps_in_files(make_settings(L, t, nr, α)) 
        @assert min_steps == N "This bundle is not fully broken! $min_steps != $N"
    end
    
    get_data(key; divide=true) = map(x -> x[key] ./(N / (divide ? 1 : N)), files_and_t)
    
    k_N = [1:n for n in N]./N
    
    seeds = round(Int64, minimum(get_data("nr_seeds_used", divide=false)))
    println("$nr, $L, $ts, $seeds")
    lables = permutedims([L"t_0 = "*"$t" for t in ts])
    colors = theme_palette(:auto)[1:length(ts)]
    
    function add_spanning_point(y_data)
        x_data = get_data("average_spanning_cluster_step")
        y = [y[round(Int64, x*N)] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, color=colors, label=nothing, markershape=:x)
    end

    nr_clusters = get_data("average_nr_clusters")
    largest_cluster = get_data("average_largest_cluster")
    largest_perimiter = get_data("average_largest_perimiter")
    most_stressed_fiber = get_data("average_most_stressed_fiber")

    yLabel(string) = use_y_lable ? string : ""

    nr_clusters_plot = plot(k_N, nr_clusters, label = lables, legend=:topright,
    ylabel=yLabel(L"\#C/N"), plot_titlevspan=-0.1, title="$nr")
    add_spanning_point(nr_clusters)

    most_stressed_fiber_plot = plot(k_N, most_stressed_fiber, label = lables, legend=:topright,
    ylabel=yLabel(L"σ"))#, title="Stress of most stressed fiber")    
    add_spanning_point(most_stressed_fiber)

    largest_cluster_plot = plot(k_N, largest_cluster, label = lables, legend=:topleft,
    ylabel=yLabel(L"S_{\mathrm{max}}/N"), plot_titlevspan=-0.1)#, title="Size of largest cluster")
    add_spanning_point(largest_cluster)

    largest_perimiter_plot = plot(k_N, largest_perimiter, label = lables, legend=:topright,
    xlabel=L"k/N", ylabel=yLabel(L"H_{\mathrm{max}}/N"))#, title="Length of the longest perimeter",)
    add_spanning_point(largest_perimiter)

    l = @layout [
        A B; C D
    ]
    plots = [nr_clusters_plot, most_stressed_fiber_plot, largest_cluster_plot, largest_perimiter_plot]
    #p = plot(plots... layout=l, plot_title="$nr "*L"L="*"$L")
    return plots
end

L = 128
ts = [0.0, 0.1, 0.2, 0.3, 0.7, 0.9]
α = 2.0
nr = ["LLS", "CLS"]
nrs = length(nr)
nr_plots = [basicPropertiesPlot(L, ts, nr[i], use_y_lable=i==1) for i in 1:nrs]
plots = reduce(vcat, reduce(vcat, collect.(zip(nr_plots...))))
plot(plots..., layout=(length(plots)÷nrs,nrs), size=(600,800), left_margin=2Plots.mm)

savefig("plots/Graphs/BundleProperties.pdf")

println("Saved plot!")
