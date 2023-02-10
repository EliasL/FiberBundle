using Plots
using JLD2
using LaTeXStrings
using Trapz
include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")

function otherPropertiesPlot(L, ts, NR; use_y_lable=true, add_ELS=true)
    
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel, labels; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, series_annotation=[], series_position=:left, log=:identity, anotate_both=true)
        if series_position==:left
            sph=:bottom
        else
            sph=:top
        end
        if !anotate_both
            series_annotation = [["" for _ in ts] text.([t in series_annotation ?  L" t_0 = "*"$t" : "" for t in ts], halign=series_position, valign=sph)]
        else
            series_annotation = permutedims([text.([t in sa ?  L" t_0 = "*"$t" : "" for t in ts], halign=series_position, valign=sph) for sa in series_annotation])
        end
        p = plot(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims,  markersize=5, xaxis=log, yaxis=log,
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, markershape=[:x :+], mc=:auto, msc=:auto,
        series_annotations =series_annotation)
        return p
    end

    function make_plot2(y, ylabel, labels; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, log=:identity)

        p = plot(x, y, label = labels, legend=false, xlims=xlims, ylims=ylims,  markersize=5, 
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=log, yaxis=log)
        return p
    end

    labels = permutedims(NR)

    largest_cluster_spanning = get_data_spanning(L, nr, ts, "average_largest_cluster", divide=:max)
    most_stressed_fiber_spanning = get_data_spanning(L, nr, ts, "average_most_stressed_fiber", divide=:max)
    println(size(largest_cluster_spanning))
    size_over_σ_LLS = make_plot(largest_cluster_spanning[:, :, 1], log=:log, series_annotation=[[],[],[0.16 0.2 0.3]], 
                        L"log$(s_{\mathrm{max}}/N)$", permutedims(["L=$l LLS" for l in L]), x=most_stressed_fiber_spanning[:, :, 1], xlabel=L"log$(σ)$", position=:topright, )
    

    size_over_σ_CLS = make_plot(largest_cluster_spanning[:, :, 2], log=:log, series_annotation=[[],[],[0.1 0.3 0.5]],
                        L"log$(s_{\mathrm{max}}/N)$", permutedims(["L=$l CLS" for l in L]), x=most_stressed_fiber_spanning[:, :, 2], xlabel=L"log$(σ)$", position=:bottomleft, series_position=:right)
    
    other_plots = [size_over_σ_LLS, size_over_σ_CLS]
    return other_plots
end

L = [32, 64, 128]
α = 2.0
nr = ["LLS", "CLS"]
ts = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts2 = (0:9) ./ 10
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr)
p = plot(plots..., size=(700,300))
savefig(p, "plots/Graphs/otherBundleProperties.svg")

println("Saved plot!")
