using Plots
using JLD2
using LaTeXStrings
using Trapz
include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")

function otherPropertiesPlot(L, ts, NR; use_y_lable=true, add_ELS=true)
    N = L.*L
    k_N = [1:n for n in N]./N
    files_nr_and_t = []
    for (i,nr) in enumerate(NR)
        push!(files_nr_and_t, [])
        for t in ts
            push!(files_nr_and_t[i], load_file(L, α, t, nr))
            # Check that the data we use is from completely borken bundles
            min_steps = get_min_steps_in_files(make_settings(L, t, nr, α)) 
            @assert min_steps == N "This bundle is not fully broken! $min_steps != $N"
        end
    end
    

    get_data(key; divide=N) = [map(x -> x[key] ./divide, f_and_t) for f_and_t in files_nr_and_t]
    get_data_x(key, x, divide=N) = [[y[round(Int64, x)] for (x,y) in zip(x[i], data_nr)] for (i,data_nr) in enumerate(get_data(key, divide=divide))]

    spanning_x = get_data("average_spanning_cluster_step", divide=1)
    get_spanning_data(key; divide=N) = get_data_x(key, spanning_x)

    most_stressed_fiber = get_data("average_most_stressed_fiber", divide=1)
    energy_x = [[argmax(σ_to_energy(σ)[2]) for σ in most_stressed_fiber_nr] for most_stressed_fiber_nr in most_stressed_fiber]
    get_energy_data(key; divide=N) = get_data_x(key, energy_x)
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel, labels; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, series_annotation=[], series_position=:left, log=:identity, anotate_both=true)
        if !anotate_both
            series_annotation = [["" for _ in ts] text.([t in series_annotation ?  L" t_0 = "*"$t" : "" for t in ts], halign=series_position, valign=:bottom)]
        else
            series_annotation = text.([t in series_annotation ?  L" t_0 = "*"$t" : "" for t in ts], halign=series_position, valign=:bottom)
        end
        p = scatter(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims,  markersize=5, xaxis=log, yaxis=log,
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

    nr_clusters = get_data("average_nr_clusters")
    largest_cluster = get_data("average_largest_cluster")
    largest_cluster_spanning = get_spanning_data("average_largest_cluster")
    largest_perimiter_spanning = get_spanning_data("average_largest_perimiter")
    most_stressed_fiber_spanning = get_spanning_data("average_most_stressed_fiber")

    perimiter_over_cluster = make_plot(largest_cluster_spanning, L"Log of spanning $s_{\mathrm{max}}/N$", labels, log=:log, anotate_both=false,
                x=largest_perimiter_spanning, xlabel=L"Log of spanning $h_{\mathrm{max}}/N$", series_annotation=[0.0, 0.2, 0.3, 0.8])
    max_σ = make_plot([maximum.(most_stressed_fiber_nr) for most_stressed_fiber_nr in most_stressed_fiber],
                        position=:topleft, L"σ_{\mathrm{max}}", labels, x=ts, xlabel=L"t_0")
#=     σ_area = make_plot([trapz.([k_N for _ in ts], most_stressed_fiber_nr) for most_stressed_fiber_nr in most_stressed_fiber],
                        position=:topleft, L"Area under $σ$", labels, x=ts, xlabel=L"t_0")
    # Scale with localization and strenght at that point
    σ_scaled = make_plot2([most_stressed_fiber_nr for most_stressed_fiber_nr in most_stressed_fiber[2]./get_energy_data("average_most_stressed_fiber", divide=1)[2]],
                        L"Scaled $σ$", permutedims(ts2), x=[k_N for _ in ts] ./ 1#= energy_x[2] =#, xlabel=L"t_0", position=:topleft,)
 =#
    size_over_σ = make_plot(largest_cluster_spanning, log=:log, series_annotation=[0.2, 0.3], series_position=:right,
                        L"log$(s_{\mathrm{max}}/N)$", labels, x=most_stressed_fiber_spanning, xlabel=L"log$(σ)$", position=:topright, )
    

    size_over_σ_LLS = make_plot(largest_cluster_spanning[1], log=:log, series_annotation=[0.1, 0.2, 0.3], 
                        L"log$(s_{\mathrm{max}}/N)$", "LLS", x=most_stressed_fiber_spanning[1], xlabel=L"log$(σ)$", position=:topright, )
    

    size_over_σ_CLS = make_plot(largest_cluster_spanning[2], log=:log, series_annotation=[0.2, 0.3, 0.4],
                        L"log$(s_{\mathrm{max}}/N)$", "CLS", x=most_stressed_fiber_spanning[2], xlabel=L"log$(σ)$", position=:topright, )
    
             
    #= σ_over_size = make_plot2(most_stressed_fiber[1], L"LLS $σ$", labels,
                        x=largest_cluster[1], xlabel=L"LLS $s_{\mathrm{max}}/N$")

    σ_over_size2 = make_plot2(most_stressed_fiber[2], L"CLS $σ$", labels,
                        x=largest_cluster[2], xlabel=L" CLS $s_{\mathrm{max}}/N$")=#

    
    
    other_plots = [size_over_σ_LLS, size_over_σ_CLS]
    return other_plots
end

L = 128
α = 2.0
nr = ["LLS", "CLS"]
ts = vcat((0:20) ./ 50, (5:9) ./ 10)
ts2 = (0:10) ./ 50
ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr)
p = plot(plots..., size=(700,300))
savefig(p, "plots/Graphs/otherBundleProperties.svg")

println("Saved plot!")
