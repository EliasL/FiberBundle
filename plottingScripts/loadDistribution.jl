using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")

function plotLoadDistribution(L, ts, nr, dist; use_y_lable=true)
    data_path="newData/"
    add_ELS=true
    N = L.*L
    files_and_t = []
    for t in ts
        push!(files_and_t, load_file(L, α, t, nr, dist, data_path=data_path))
        # Check that the data we use is from completely borken bundles
        min_steps = get_min_steps_in_files(make_settings(L, t, nr, α, dist, data_path)) 
        @assert min_steps == N "This bundle, $nr L=$L t0=$t, is not fully broken! $min_steps != $N"
    end

    if nr=="LLS" && add_ELS
        push!(files_and_t, load_file(L, α, 0.5, "ELS", dist, data_path=data_path))
        # Check that the data we use is from completely borken bundles
        min_steps = get_min_steps_in_files(make_settings(L, 0.0, "ELS", α, dist, data_path)) 
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
    colors = vcat(theme_palette(:auto)[1:12],theme_palette(:auto)[1:12])[1:length(labels)]
    
    function add_points(y_data)
        # Add spanning point
        x_data = get_data("average_spanning_cluster_step")
        y = [y[round(Int64, x*N)] for (x,y) in zip(x_data,y_data)]
        #Draw spanning
        scatter!(x_data, y, color=colors, label=nothing, markershape=:x, markeralpha=1) 

       #=  #Add energy change point
        x_data = [argmax(σ_to_energy(σ)[2]) for σ in most_stressed_fiber]    
        y = [y[round(Int64, x)] for (x,y) in zip(x_data,y_data)]
        #Draw localization
        scatter!(x_data/N, y, color=colors, label=nothing, markershape=:vline, markersize=10, markeralpha=1, markerstrokewidth=1)
        =# 
        #Add max σ_max
        x_data = [argmax(σ) for σ in most_stressed_fiber]    
        y = [y[round(Int64, x)] for (x,y) in zip(x_data,y_data)]
        #Draw localization
        scatter!(x_data/N, y, markerstrokecolor=colors, markercolor=:transparent, label=nothing, markershape=:diamond, markersize=5, markerstrokewidth=1)

    end
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(0, 1.3), position=:topright)
        # Use empty scatter as title
        plot = scatter([0],[0], label=L"t_0", ms=0, mc=:white, msc=:white)
        plot!(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims, color= permutedims(colors),
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, yaxis=:identity,
        linestyle=hcat([:dash], permutedims([:solid for _ in 1:(length(ts)-(add_ELS ? 1 : 2))]),[:dot]))
        add_points(y)
        return plot
    end

    nr_clusters = get_data("average_nr_clusters")
    largest_cluster = get_data("average_largest_cluster")
    largest_perimiter = get_data("average_largest_perimiter")
    most_stressed_fiber = get_data("average_most_stressed_fiber", divide=1)

    load = []
    for i in eachindex(largest_cluster)
        push!(load, [])
        for j in eachindex(largest_cluster[i])
            v = largest_perimiter[i][j] / largest_cluster[i][j] #/ largest_perimiter[i][j] * nr_clusters[i][j]
            push!(load[i], v)
        end
    end

    return make_plot(load, L"\frac{\langle s \rangle}{\langle h \rangle} M", ylims=(-Inf, 0.23))
end

L = 128
ts = round.((1 .- vcat((0:20) ./ 50, (5:7) ./ 10)) ./2, digits=2)
ts = vcat(0.05:0.05:0.25, 0.3:0.01:0.5)
ts = [0.1, 0.27, 0.3, 0.35, 0.40, 0.5]
α = 2.0
nr = ["LLS", "CLS"]
dist = "ConstantAverageUniform"
nrs = length(nr)
plots = [plotLoadDistribution(L, ts, nr[i], dist, use_y_lable=i==1) for i in 1:nrs]
p = plot(plots..., layout=(1,nrs), size=(300,200), left_margin=2Plots.mm)

savefig(p, "plots/Graphs/$(dist)_LoadDist.svg")

println("Saved plot!")
