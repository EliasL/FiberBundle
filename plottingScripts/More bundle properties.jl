using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")

function otherPropertiesPlot(L, ts, nr; use_y_lable=true, basic=true, add_ELS=true)
    N = L.*L
    files_and_t = []
    for t in ts
        push!(files_and_t, load_file(L, α, t, nr))
        # Check that the data we use is from completely borken bundles
        min_steps = get_min_steps_in_files(make_settings(L, t, nr, α)) 
        @assert min_steps == N "This bundle is not fully broken! $min_steps != $N"
    end
    
    get_data(key; divide=N) = map(x -> x[key] ./divide, files_and_t)

    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(0, 1.2), possition=:topright)
        # Use empty scatter as title
        if basic
            plot = scatter([0],[0], label=L"t_0", ms=0, mc=:white, msc=:white)
        else
            plot = scatter([0],[0], label="", ms=0, mc=:white, msc=:white)
        end
        if basic
            plot!(x, y, label = labels, legend=possition, xlims=xlims, ylims=ylims, color= permutedims(colors),
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title,
            linestyle=hcat([:dash], permutedims([:solid for _ in 1:(length(ts)-1)]),[:dot]))
        else
            plot!(x, y, label = labels, legend=possition,xlabel=xlabel, ylabel=yLabel(ylabel), title=title,)
        end
        if x == k_N
            add_points(y)
        end
        return plot
    end

    nr_clusters = get_data("average_nr_clusters")
    largest_cluster = get_data("average_largest_cluster")
    largest_perimiter = get_data("average_largest_perimiter")
    most_stressed_fiber = get_data("average_most_stressed_fiber", divide=1)


    if basic
        nr_clusters_plot = make_plot(nr_clusters, L"M/N", title=nr*(nr=="LLS" && add_ELS ? " and ELS" : ""), ylims=(0,0.13))
        most_stressed_fiber_plot = make_plot(most_stressed_fiber,L"σ")
        largest_cluster_plot = make_plot(largest_cluster,L"s_{\mathrm{max}}/N")    
        largest_perimiter_plot = make_plot(largest_perimiter,L"h_{\mathrm{max}}/N", ylims=(0,0.375), xlabel=L"k/N")
        basic_plots = [nr_clusters_plot, most_stressed_fiber_plot, largest_cluster_plot, largest_perimiter_plot]
        return basic_plots
    else 

        spanning_x = get_data("average_spanning_cluster_step")
        largest_cluster_spanning = [y[round(Int64, x*N)] for (x,y) in zip(spanning_x,largest_cluster)]
        largest_perimiter_spanning= [y[round(Int64, x*N)] for (x,y) in zip(spanning_x,largest_perimiter)]
        perimiter_over_cluster = make_plot(largest_cluster_spanning, L"Spanning $s_{\mathrm{max}}/N$", x=largest_perimiter_spanning, xlabel=L"Spanning $h_{\mathrm{max}}/N$")
        max_σ = make_plot(maximum.(most_stressed_fiber), L"σ_{\mathrm{max}}", x=ts, xlabel=L"t_0")
        other_plots = [perimiter_over_cluster, max_σ]
        return other_plots
    end
end

L = 128
α = 2.0
nr = ["LLS", "CLS"]
ts = vcat((0:20) ./ 50, (5:9) ./ 10)
plots = otherPropertiesPlot(L, ts, nr)
p = plot(plots..., size=(700,800))
savefig(p, "plots/Graphs/otherBundleProperties.svg")

println("Saved plot!")
