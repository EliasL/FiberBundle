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

    if nr=="LLS"
        push!(files_and_t, load_file(L, α, 0.0, "ELS"))
        # Check that the data we use is from completely borken bundles
        min_steps = get_min_steps_in_files(make_settings(L, 0.0, "ELS", α)) 
        @assert min_steps == N "This bundle is not fully broken! $min_steps != $N"
    end

    
    get_data(key; divide=N) = map(x -> x[key] ./divide, files_and_t)
    
    k_N = [1:n for n in N]./N
    
    seeds = round(Int64, minimum(get_data("nr_seeds_used", divide=1)))
    if nr=="LLS"
        extra_label=["ELS"]
    else
        extra_label=[]
    end
    labels = permutedims(vcat(["$t" for t in ts], extra_label))
    colors = theme_palette(:auto)[1:length(labels)]
    
    function add_spanning_point(y_data)
        x_data = get_data("average_spanning_cluster_step")
        y = [y[round(Int64, x*N)] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, color=colors, label=nothing, markershape=:x)
    end

    function make_plot(y, possition, ylabel,  xlims=(-Inf, Inf), title="", xlabel="")
        # Use empty scatter as title
        plot = scatter([0],[0], label=L"t_0", ms=0, mc=:white, msc=:white)
        plot!(k_N, y, label = labels, legend=possition, xlims=xlims,color= permutedims(colors),
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, linestyle=hcat([:dash], permutedims([:solid for _ in 1:(length(ts)-1)]),[:dot]))
        add_spanning_point(y)
        return plot
    end

    nr_clusters = get_data("average_nr_clusters")
    largest_cluster = get_data("average_largest_cluster")
    largest_perimiter = get_data("average_largest_perimiter")
    most_stressed_fiber = get_data("average_most_stressed_fiber", divide=1)

    yLabel(string) = use_y_lable ? string : ""

    nr_clusters_plot = make_plot(nr_clusters, :bottomright, L"\#C/N", (0,1.249), nr*(nr=="LLS" ? "/ELS" : ""))

    most_stressed_fiber_plot = make_plot(most_stressed_fiber, :bottomright,L"σ_{\mathrm{max}}", (0,1.249),)

    largest_cluster_plot = make_plot(largest_cluster, :bottomright,L"S_{\mathrm{max}}/N", (0,1.249),)
    
    largest_perimiter_plot = make_plot(largest_perimiter, :bottomright,L"H_{\mathrm{max}}/N", (0,1.249), "", L"k/N")
    
    l = @layout [
        A B; C D
    ]
    plots = [nr_clusters_plot, most_stressed_fiber_plot, largest_cluster_plot, largest_perimiter_plot]
    #p = plot(plots... layout=l, plot_title="$nr "*L"L="*"$L")
    return plots
end

L = 128
ts = [0.0, 0.1, 0.2, 0.3, 0.7]
α = 2.0
nr = ["LLS", "CLS"]
nrs = length(nr)
nr_plots = [basicPropertiesPlot(L, ts, nr[i], use_y_lable=i==1) for i in 1:nrs]
plots = reduce(vcat, reduce(vcat, collect.(zip(nr_plots...))))
plot(plots..., layout=(length(plots)÷nrs,nrs), size=(700,800), left_margin=2Plots.mm)

savefig("plots/Graphs/BundleProperties.pdf")

println("Saved plot!")
