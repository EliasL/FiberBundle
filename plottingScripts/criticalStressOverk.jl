using Plots
using JLD2
using LaTeXStrings
using Trapz
using LsqFit
include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")


function otherPropertiesPlot(L, ts, NR, dist; use_y_lable=true)
    yLabel(string) = use_y_lable ? string : ""

    colors = permutedims(theme_palette(:auto)[1:16])
    markershape = [:+ :star4 :diamond :rect :ltriangle :rtriangle :star6]

    function make_plot1(X, Y, ylabel, nr_index; 
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf), label=permutedims(L),
        position=:topright, log_scale=:identity, title="")
        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=:identity,
            yaxis=log_scale, framestyle=:box)
        #plot!([], [], label=L"L", alpha=0)

        clusterSize = get_data_kN(L, nr, ts, dist, "average_largest_cluster", return_kN=false)
        localization = [find_localization(clusterSize[1][:, i, nr_index], critical_gradient)/length(clusterSize[1][:, i, 1]) for i in eachindex(ts)]
        add_plot(ts, localization, 2, label=L"k_l")
        
        localization = [find_localization(clusterSize[1][:, i, nr_index], 0.5)/length(clusterSize[1][:, i, 1]) for i in eachindex(ts)]
        add_plot_not_scatter(ts, localization, 3, label=L"k_l \: (1/2)")
    
        scatter!(X, Y, label=label,markeralpha=1.0, c=colors[1],
        legend=position, markershape=markershape, markersize=7,
        markerstrokecolor=colors[1], z_order=:back)
        
        if dist=="Weibull"
            xflip!(p)
        end
        return p
    end

    function add_plot(X, Y, i; label="")
        scatter!(X, Y, label=label, linestyle=:solid, markersize=7,
        markerstrokecolor=colors[i],
        markershape=markershape[i], z_order=:front)
    end
    function add_plot_not_scatter(X, Y, i; label="")
        plot!(X, Y, label=label, linestyle=:solid, markersize=7, c=colors[i],
        markerstrokecolor=colors[i], z_order=:front)
    end
    labels = permutedims(NR)

    σ_c, x = get_data(L, nr, ts, dist, "most_stressed_fiber",
        "most_stressed_fiber", argmax, ex=[0, 0], rel_x=true, average=false, return_x=true,
        data_path=data_path)
    
        
    LLS_σ_c_N_plot = make_plot1(ts, x[:, :, 1], label=L"k_c",
        L"k/N", 1, log_scale=:identity, title="LLS",
        xlabel=xlabel, position=pos,)

#=     most_stressed_fiber = get_data_kN(L, nr, ts, dist, "average_most_stressed_fiber", return_kN=false)
    E = [σ_to_energy(most_stressed_fiber[1][:, i, 1]) for i in eachindex(ts)]
    k = [argmax(E[1])/length(E[1]) for E in E]
    add_plot(ts, k, label="Max energy") =#

    
    #localization = [find_localization(clusterSize[1][:, i, 1], 1/2)/length(clusterSize[1][:, i, 1]) for i in eachindex(ts)]
    #add_plot(ts, localization, label=L"k_l")
    
    CLS_σ_c_N_plot = make_plot1(ts, x[:, :, 2], label=L"k_c",
        L"k/N", 2, log_scale=:identity, title="CLS", 
        xlabel=xlabel, position=pos,)

    return [LLS_σ_c_N_plot,CLS_σ_c_N_plot]
end

L = [128]
α = 2.0
nr = ["LLS", "CLS"]

#ts = vcat((0:20) ./ 50, (5:9) ./ 10)
dist = "ConstantAverageUniform"
dist = "Weibull"
if dist=="Weibull"
    xlabel=L"t_w"
    ts = vcat([0.5], 1:0.1:1.5, 2:0.5:5)
else
    xlabel=L"t_0"
    ts = vcat(0.05:0.05:0.20, 0.25:0.01:0.5)
end

#ts = (0:7) ./ 10
data_path = "newData/"
critical_gradient = 1
pos = :topleft
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr, dist)
xpsize = 230
ypsize = 230
println("plot")
p = plot(plots..., size=(xpsize * 1.1 * length(plots), ypsize *
                                                       maximum([length(plots) / 2, 1])))

#p2 = plot(plots[3:4]..., size=(psize*length(nr)*1.1,psize
#length(plots)/2/length(nr)), layout = @layout([ A B;]))
display(p)
savefig(p, "plots/Graphs/$(dist)_CriticalStressOverk.pdf")
#savefig(p2, "plots/Graphs/$(dist)_s_over_sigma.pdf")

println("Saved plot!")

