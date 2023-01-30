using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")

function basicPropertiesPlot(L, ts, nr; use_y_lable=true)
    add_ELS=true
    N = L.*L
    files_and_t = []
    for t in ts
        push!(files_and_t, load_file(L, α, t, nr))
        # Check that the data we use is from completely borken bundles
        min_steps = get_min_steps_in_files(make_settings(L, t, nr, α)) 
        @assert min_steps == N "This bundle is not fully broken! $min_steps != $N"
    end

    if nr=="LLS" && add_ELS
        push!(files_and_t, load_file(L, α, 0.0, "ELS"))
        # Check that the data we use is from completely borken bundles
        min_steps = get_min_steps_in_files(make_settings(L, 0.0, "ELS", α)) 
        @assert min_steps == N "This bundle is not fully broken! $min_steps != $N"
    end

    
    get_data(key; divide=N) = map(x -> x[key] ./divide, files_and_t)
    
    k_N = [1:n for n in N]./N
    
    seeds = round(Int64, minimum(get_data("nr_seeds_used", divide=1)))
    if nr=="LLS" && add_ELS
        extra_label=["ELS"]
    else
        extra_label=[]
    end
    labels = permutedims(vcat(["$t" for t in ts], extra_label))
    colors = theme_palette(:auto)[1:length(labels)]
    
    function add_points(y_data)
        # Add spanning point
        x_data = get_data("average_spanning_cluster_step")
        y = [y[round(Int64, x*N)] for (x,y) in zip(x_data,y_data)]
        #Draw spanning
        scatter!(x_data, y, color=colors, label=nothing, markershape=:x) 

        #Add energy change point
        x_data = [argmax(σ_to_energy(σ)[1]) for σ in most_stressed_fiber]    
        y = [y[round(Int64, x)] for (x,y) in zip(x_data,y_data)]
        #Draw localization
        scatter!(x_data/N, y, color=colors, label=nothing, markershape=:vline, markersize=10)
        
        #Add max σ_max
        x_data = [argmax(σ) for σ in most_stressed_fiber]    
        println(x_data)
        y = [y[round(Int64, x)] for (x,y) in zip(x_data,y_data)]
        #Draw localization
        scatter!(x_data/N, y, markerstrokecolor=colors, markercolor=:transparent, label=nothing, markershape=:diamond, markersize=5, markerstrokewidth=1)


        #= # Add localization point
        s_data = get_data("average_largest_cluster")
        function large_slope(s)
            return s>1/N
        end
        function large_cluster(s,t)
            return s.>0.0000002*N*(1-t)
        end
        x_data = [ findfirst(large_slope, diff(s))*(1-t) for (s,t) in zip(s_data,ts)]
        #x_data = [ findfirst(large_cluster(s,t)) for (s,t) in zip(s_data,ts)]
        y = [y[round(Int64, x)] for (x,y) in zip(x_data,y_data)]
        #Draw localization
        scatter!(x_data/N, y, color=colors, label=nothing, markershape=:+) =#
    end
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(0, 1.2), position=:topright)
        # Use empty scatter as title
        plot = scatter([0],[0], label=L"t_0", ms=0, mc=:white, msc=:white)
        plot!(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims, color= permutedims(colors),
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title,
        linestyle=hcat([:dash], permutedims([:solid for _ in 1:(length(ts)-1)]),[:dot]))
        add_points(y)
        return plot
    end

    nr_clusters = get_data("average_nr_clusters")
    largest_cluster = get_data("average_largest_cluster")
    largest_perimiter = get_data("average_largest_perimiter")
    most_stressed_fiber = get_data("average_most_stressed_fiber", divide=1)


    nr_clusters_plot = make_plot(nr_clusters, L"M/N", title=nr*(nr=="LLS" && add_ELS ? " and ELS" : ""), ylims=(0,0.13))
    most_stressed_fiber_plot = make_plot(most_stressed_fiber,L"σ")
    largest_cluster_plot = make_plot(largest_cluster,L"s_{\mathrm{max}}/N")    
    largest_perimiter_plot = make_plot(largest_perimiter,L"h_{\mathrm{max}}/N", ylims=(0,0.375), xlabel=L"k/N")
    basic_plots = [nr_clusters_plot, most_stressed_fiber_plot, largest_cluster_plot, largest_perimiter_plot]
    return basic_plots
end

L = 128
ts = [0.0, 0.1, 0.2, 0.3, 0.7]
α = 2.0
nr = ["LLS", "CLS"]
nrs = length(nr)
nr_plots = [basicPropertiesPlot(L, ts, nr[i], use_y_lable=i==1) for i in 1:nrs]
plots = reduce(vcat, reduce(vcat, collect.(zip(nr_plots...))))
p = plot(plots..., layout=(length(plots)÷nrs,nrs), size=(700,800), left_margin=2Plots.mm)

savefig(p, "plots/Graphs/BundleProperties.svg")

println("Saved plot!")
