using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")


function maxify(σ, nr_seeds, max_range, L)
    max_σ = zeros(Float64, L^2)    
    for seed in 1:nr_seeds
        current_max = 0
        current_max_k = 0
        for k in 1:L^2 
            new_value = σ[k,seed]
            if current_max < new_value
                current_max = new_value
                current_max_k = k
            end
            if k-current_max_k > max_range/2
                search_interval = k:minimum([k+round(Int64,max_range/2), L^2])
                current_max_k = argmax((view(σ,search_interval, seed)))
                current_max = σ[search_interval[current_max_k], seed]
            end
            max_σ[k] += current_max
        end
    end 
    max_σ ./= nr_seeds     
    return max_σ
end


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

    
    simple_get_data(key; divide=N) = map(x -> x[key] ./divide, files_and_t)
    
    k_N = [1:n for n in N]./N
    
    seeds = round(Int64, minimum(simple_get_data("nr_seeds_used", divide=1)))
    if nr=="LLS" && add_ELS
        extra_label=["ELS"]
    else
        extra_label=[]
    end
    labels = permutedims(vcat(["$t" for t in ts], extra_label))
    colors = vcat(theme_palette(:auto)[1:12],theme_palette(:auto)[1:12])[1:length(labels)]
    
    function add_points(y_data)
        # Draw spanning
#=         x_data = simple_get_data("average_spanning_cluster_step")
        println(size(x_data))
        y = [y[round(Int64, x*length(y))] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, color=colors, label=nothing,
                markersize=6, markershape=:vline, markeralpha=1)  =#

        # Draw average max stress
#=         x_data = [argmax(σ)/length(σ) for σ in most_stressed_fiber]   
        y = [y[round(Int64, x*length(y))] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, markerstrokecolor=colors, markercolor=:transparent,
        label=nothing, markershape=:diamond, markersize=5, markerstrokewidth=1) =#

        # Draw critical stressed
        x_data = σ_c_x[:, 1, 1]
        y = [y[round(Int64, x*length(y))] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, markerstrokecolor=colors, markercolor=:transparent,
        label=nothing, markershape=:diamond, markersize=6, markerstrokewidth=1)


        # Draw localization
        
        clusterSize = get_data_kN(L, [nr], ts, dist, "average_largest_cluster", return_kN=false)
        x_data = [find_localization(clusterSize[1][:, i, 1])/length(clusterSize[1][:, i, 1]) for i in eachindex(ts)]
        y = [y[round(Int64, x*length(y))] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, c=colors, markeralpha=1, 
        label=nothing, markershape=:+, markersize=8, markerstrokewidth=1)

        #= nrClusters = get_data_kN(L, [nr], ts, dist, "average_nr_clusters", return_kN=false)
        x_data = [find_localization_nr_clusters(nrClusters[1][:, i, 1])/length(nrClusters[1][:, i, 1]) for i in eachindex(ts)]
        y = [y[round(Int64, x*length(y))] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, c=colors, markeralpha=1, 
        label=nothing, markershape=:+, markersize=5, markerstrokewidth=1) =#

#=         threshhold=0.01
        x_data = [findfirst(x->x>threshhold, s)/length(s) for s in largestluster]
        y = [y[round(Int64, x*length(y))] for (x,y) in zip(x_data,y_data)]
        scatter!(x_data, y, c=colors, markeralpha=1, 
        label=nothing, markershape=:vline, markersize=8, markerstrokewidth=1) =#
        
    end
    


        function get_σ_data(ts)
            Y = []
            max_range = L^2/800
            nr_seeds = 1000
            for (i, t) in zip(1:length(ts), ts)
                println(t)
                σ, x = get_data_kN(L, [nr], t, dist, "most_stressed_fiber",
                average=false, return_kN=true, divide=1, nr_seeds=nr_seeds)
                σ = σ[1][:, 1, 1, :]
                
                #push!(Y, maxify(σ, nr_seeds, max_range, L))
                push!(Y,vec(mean(σ, dims=2)))
            end
            return Y
        end
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(0, 1.35), position=:topright)
        # Use empty scatter as title
        #plot = scatter([0.01], label=L"t_0", ms=0, mc=:white, msc=:white)
        p=plot()
        label = dist=="Weibull" ? L"t_w" : L"t_0"
        plot!([0.01], label=label, ms=0, mc=:white, msc=:white, c=:white)
        plot!(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims, color= permutedims(colors),
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, size=(300,250),
        linestyle=hcat([:dash], permutedims([:solid for _ in 1:(length(ts)-(add_ELS ? 1 : 2))]),[:dot]))
        add_points(y)
        return p
    end

    nrlusters = simple_get_data("average_nr_clusters")
    largestluster = simple_get_data("average_largest_cluster")
    largest_perimiter = simple_get_data("average_largest_perimiter")
    σ_c, σ_c_x = get_data(L, [nr], ts, dist, "most_stressed_fiber",
        "most_stressed_fiber", argmax, ex=[0, 0], rel_x=true, average=false,
        return_x=true, data_path=data_path)
    most_stressed_fiber = get_σ_data(ts)
    
#=     println(size(most_stressed_fiber[1]))
    println(size(largest_perimiter[1])) =#
    #most_stressed_fiber = get_data("average_most_stressed_fiber", divide=1)
    most_stressed_fiber_plot = make_plot(most_stressed_fiber,L"\langle σ \rangle", title=nr*(nr=="LLS" && add_ELS ? " and ELS" : ""))
    #most_stressed_fiber_plot = make_none_averaged_σ_plot(L, ts, nr, dist)
    largestluster_plot = make_plot(largestluster,L"\langle s_{\mathrm{max}}/N \rangle")    
    largest_perimiter_plot = make_plot(largest_perimiter,L"\langle h_{\mathrm{max}}/N \rangle", ylims=(0,0.375))
    nrlusters_plot = make_plot(nrlusters, L"\langle M/N \rangle", ylims=(0,0.13), xlabel=L"k/N")
    basic_plots = [most_stressed_fiber_plot, largestluster_plot, largest_perimiter_plot, nrlusters_plot]
    return basic_plots
end



function make_none_averaged_σ_plot(L, ts, nr, dist, pos)
    nr_seeds = 100
    divisions=L*L
    p = plot()
    colors = theme_palette(:auto)[1:length(ts)]
    for (i, t, c) in zip(1:length(ts), ts, colors)
        raw_σ, x = get_data_kN(L, nr, t, dist, "most_stressed_fiber",
        average=false, return_kN=true, divide=1, nr_seeds=nr_seeds)
        every=round(Int64, L^2/divisions)
        
        σ = collect([map(maximum, Iterators.partition(raw_σ[1][:, 1, 1, i], every)) for i in 1:nr_seeds])
        #σ = mean(σ, dims=1)
        x = collect([map(maximum, Iterators.partition(x[1], every)) for i in 1:nr_seeds])
        
        p = plot!(x,  σ, label="", title=nr[1],
        xlabel=L"k/N", ylabel=(nr[1]=="CLS" ? "" : "max("*L"σ)_{100}"), c=c, alpha=0.05)
        #add_points(mean(raw_σ[1][:, 1, 1, :], dims=2), i=i)
    end
    plot!([1.0], label=L"t_0", ms=0, mc=:white, msc=:white, xlims=xlims(p), ylims=ylims(p), c=:white)
    plot!(ones((1, length(ts))), label=permutedims(ts), c=permutedims(colors), xlims=xlims(p), ylims=ylims(p), legend=pos)
    return p
end

L = 128
ts = round.((1 .- vcat((0:20) ./ 50, (5:7) ./ 10)) ./2, digits=2)
ts = vcat(0.05:0.05:0.25, 0.3:0.01:0.5)
α = 2.0
nr = ["LLS", "CLS"]
dist = "Weibull"
dist = "ConstantAverageUniform"
if dist == "Weibull"
    ts = reverse([0.5, 1, 1.5, 2, 3, 4, 5])
else
    ts = [0.1, 0.27, 0.3, 0.35, 0.40, 0.5]
end
#= nrs = length(nr)
println("Started...")
nr_plots = [basicPropertiesPlot(L, ts, nr[i], dist, use_y_lable=i==1) for i in 1:nrs]
plots = reduce(vcat, reduce(vcat, collect.(zip(nr_plots...))))
names = ["sigma", "cluster_size", "perimiter_length", "nrClusters"]
yValues = [L"\langle σ \rangle", L"\langle s_{\mathrm{max}}/N \rangle",L"\langle h_{\mathrm{max}}/N \rangle", L"\langle M/N \rangle"] 
for i in eachindex(plots)
    p = plots[i]
    p = plot(p,  xlabel=L"k/N", ylabel=yValues[ceil(Int64,i/2)], size=(340, 220), title="")
    savefig(p, "plots/Graphs/Basic/$(dist) $(nr[mod1(i, 2)]) $(names[ceil(Int64,i/2)]).pdf")
end
p = plot(plots..., layout=(length(plots)÷nrs,nrs), size=(700,800), left_margin=2Plots.mm, link=:x)
savefig(p, "plots/Graphs/$(dist)_BundleProperties.pdf") =#

L=128
nr = ["LLS", "CLS"]
ts = [0.5, 0.4, 0.3, 0.27]
Plots.scalefontsizes(4)
p = [make_none_averaged_σ_plot(L, ts, [nr], dist, (nr=="LLS" ? :bottom : :topright)) for nr=nr]
#= p1 = make_none_averaged_σ_plot(L, [0.4], ["LLS"], dist)
p2 = make_none_averaged_σ_plot(L, [0.4], ["CLS"], dist) =#
p = plot(p..., size=(1200*length(nr), 1200), layout=(1,length(nr)))
savefig(p, "plots/Graphs/$(dist)_non_averaged_NotMAX.png")

dist = "Weibull"
ts = [0.5, 1, 1.5, 2]

p = [make_none_averaged_σ_plot(L, ts, [nr], dist, (nr=="LLS" ? :bottom : :topright)) for nr=nr]
#= p1 = make_none_averaged_σ_plot(L, [0.4], ["LLS"], dist)
p2 = make_none_averaged_σ_plot(L, [0.4], ["CLS"], dist) =#
p = plot(p..., size=(1200*length(nr), 1200), layout=(1,length(nr)))
savefig(p, "plots/Graphs/$(dist)_non_averaged_NotMAX.png")
println("Saved plot!")
