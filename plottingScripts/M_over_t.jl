using Plots
using JLD2
using LaTeXStrings
using Trapz
using LsqFit

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")


function otherPropertiesPlot(L, ts, NR, dist; use_y_lable=true, add_ELS=true)
    colors = permutedims(theme_palette(:auto)[1:16])
    markershape = [:circle :star4 :diamond :rect :ltriangle :rtriangle :star6]


    function make_plot1(X, Y, ylabel; 
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf), label=permutedims(L),
        position=:topleft, log_scale=:identity, title="")
        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=ylabel, title=title, xaxis=:identity,
            yaxis=log_scale, framestyle=:box)
        #plot!([], [], label=L"L", alpha=0)
        plot!([0.01], [0], label="L", ms=0, mc=:white, msc=:white, c=:white)
        scatter!(X, Y, label=label, linestyle=:solid,
        legend=position, markershape=markershape, markersize=7,
        markerstrokecolor=colors, framestyle="")
        extra_plot_nr=2
        return p
    end

    labels = permutedims(NR)
    nr_clusters = get_data_kN(L, NR, ts, dist, "average_nr_clusters", return_kN=false)
    
    max_nr_clusters = maximum.(nr_clusters, dims=1)
    max_nr_clust_LLS = zeros(length(ts), length(L))
    max_nr_clust_CLS = zeros(length(ts), length(L))
    for t=eachindex(ts), l=eachindex(L) 
      max_nr_clust_LLS[t, l] = max_nr_clusters[l][1, t, 1] 
    end
    for t=eachindex(ts), l=eachindex(L) 
      max_nr_clust_CLS[t, l] = max_nr_clusters[l][1, t, 2] 
    end
    
    
    max_nr_clusters_LLS_plot = make_plot1(repeat(ts,1,length(L)), max_nr_clust_LLS,
        L"M", title="LLS",xlabel=L"t_0")
    max_nr_clusters_CLS_plot = make_plot1(repeat(ts,1,length(L)), max_nr_clust_CLS,
        "", title="CLS",xlabel=L"t_0")
    
    return [max_nr_clusters_LLS_plot,max_nr_clusters_CLS_plot ]
end

L = [32, 64, 128, 256, 512]
α = 2.0
nr = ["LLS", "CLS"]
dist = "ConstantAverageUniform"
data_path="newData/"
ts = vcat(0.05:0.05:0.25, 0.3:0.01:0.5)
#ts = vcat(0.3:0.01:0.5)
#ts2 = (0:9) ./ 10
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr, dist)
p = plot(plots..., size=(440,220))
savefig(p, "plots/Graphs/$(dist)_M_over_t0.pdf")

println("Saved plot!")
p
