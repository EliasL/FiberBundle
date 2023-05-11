using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")

function basicPropertiesPlot(L, ts, nr, dist; use_y_lable=true)
    data_path="newData/"
    add_ELS=false
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
        println(size(y_data))
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


        #= # Add localization point
        s_data = get_data("average_largestluster")
        function large_slope(s)
            return s>1/N
        end
        function largeluster(s,t)
            return s.>0.0000002*N*(1-t)
        end
        x_data = [ findfirst(large_slope, diff(s))*(1-t) for (s,t) in zip(s_data,ts)]
        #x_data = [ findfirst(largeluster(s,t)) for (s,t) in zip(s_data,ts)]
        y = [y[round(Int64, x)] for (x,y) in zip(x_data,y_data)]
        #Draw localization
        scatter!(x_data/N, y, color=colors, label=nothing, markershape=:+) =#
    end
    
    function make_none_averaged_σ_plot(L, ts, nr, dist)
            nr_seeds = 3000
            divisions=200
            p = plot()
            colors = theme_palette(:auto)[1:length(ts)]
            for (t, c) in zip(ts, colors)
                raw_σ, x = get_data_kN(L, [nr], t, dist, "most_stressed_fiber",
                average=false, return_kN=true, divide=1, nr_seeds=nr_seeds)
                every=round(Int64, L^2/divisions)
                
                σ = collect([map(maximum, Iterators.partition(raw_σ[1][:, 1, 1, i], every)) for i in 1:nr_seeds])
                σ = mean(σ, dims=1)
                x = collect([map(maximum, Iterators.partition(x[1], every)) for i in 1:nr_seeds])
                
                p = plot!(x,  σ, label="", title=nr[1]*" "*L"t_0 = "*"$(ts[1])",
                    xlabel=L"k/N", ylabel="σ", c=c, alpha=0.05)
                add_points(mean(raw_σ[1][:, 1, 1, :], dims=2)[1])
            end
            return p
    end
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(0, 1.3), position=:topright)
        # Use empty scatter as title
        plot = scatter([0],[0], label=L"t_0", ms=0, mc=:white, msc=:white)
        plot!(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims, color= permutedims(colors),
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title,
        linestyle=hcat([:dash], permutedims([:solid for _ in 1:(length(ts)-(add_ELS ? 1 : 2))]),[:dot]))
        add_points(y)
        return plot
    end

    nrlusters = get_data("average_nr_clusters")
    largestluster = get_data("average_largest_cluster")
    largest_perimiter = get_data("average_largest_perimiter")
    most_stressed_fiber = get_data("average_most_stressed_fiber", divide=1)


#    most_stressed_fiber_plot = make_plot(most_stressed_fiber,L"\langle σ \rangle", title=nr*(nr=="LLS" && add_ELS ? " and ELS" : ""))
    most_stressed_fiber_plot = make_none_averaged_σ_plot(L, ts, nr, dist)
    largestluster_plot = make_plot(largestluster,L"\langle s_{\mathrm{max}}/N \rangle")    
    largest_perimiter_plot = make_plot(largest_perimiter,L"\langle h_{\mathrm{max}}/N \rangle", ylims=(0,0.375), xlabel=L"k/N")
    nrlusters_plot = make_plot(nrlusters, L"\langle M/N \rangle", ylims=(0,0.13))
    basic_plots = [most_stressed_fiber_plot, largestluster_plot, largest_perimiter_plot, nrlusters_plot]
    return basic_plots
end


L = 128
ts = round.((1 .- vcat((0:20) ./ 50, (5:7) ./ 10)) ./2, digits=2)
ts = vcat(0.05:0.05:0.25, 0.3:0.01:0.5)
ts = [0.1, 0.27, 0.3, 0.35, 0.40, 0.5]
α = 2.0
nr = ["LLS", "CLS"]
dist = "ConstantAverageUniform"
nrs = length(nr)
nr_plots = [basicPropertiesPlot(L, ts, nr[i], dist, use_y_lable=i==1) for i in 1:nrs]
plots = reduce(vcat, reduce(vcat, collect.(zip(nr_plots...))))
p = plot(plots..., layout=(length(plots)÷nrs,nrs), size=(700,800), left_margin=2Plots.mm, link=:x)

savefig(p, "plots/Graphs/$(dist)_BundleProperties.svg")
#= L=512
ts = [0.5, 0.4, 0.3]
p = [make_none_averaged_σ_plot(L, ts, [nr], dist) for nr=nr]
#= p1 = make_none_averaged_σ_plot(L, [0.4], ["LLS"], dist)
p2 = make_none_averaged_σ_plot(L, [0.4], ["CLS"], dist) =#
p = plot(p..., size=(300*length(nr), 300), layout=(1,length(nr)))
savefig(p, "plots/Graphs/non_averaged.pdf") =#
println("Saved plot!")
