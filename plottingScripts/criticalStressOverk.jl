using Plots
using JLD2
using LaTeXStrings
using Trapz
using LsqFit
include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")

funk(X, p) = p[2] .+ p[1] .* X

function get_fit(x, y)
    p0 = [0.22, 0.22]
    return curve_fit(funk, x, y, p0)
end

function lin(x)
    return minimum(x):0.1:maximum(x)
end

function myerror(x, y, f)
    fit = funk(x, f.param)
    return sum(abs.(fit.-y))./length(y)
end

function myfit(x, y; fit_interval=1)
    r = round(Int64, fit_interval * length(x))
    f = get_fit(x[1:r], y[1:r])
    return funk(lin(x), f.param), f.param, myerror(x, y, f)
end

function otherPropertiesPlot(L, ts, NR, dist; use_y_lable=true)
    yLabel(string) = use_y_lable ? string : ""

    colors = permutedims(theme_palette(:auto)[1:16])
    markershape = [:circle :star4 :diamond :rect :ltriangle :rtriangle :star6]

    function add_fits(x, y, plot_index)
        params = zeros(Float64, (size(y)[1],2))
        errors = zeros(Float64, size(y)[1])
        nr_color = 1
        for i in 1:size(y)[1] 
            fit, param, error = myfit(x, y[i, :])
            if i in plot_index
                plot!(lin(x), fit, label="", c=colors[nr_color],
                markershape=markershape, markerstrokecolor=colors[nr_color],
                markerstrokealpha=0)
                nr_color += 1
            end
            params[i, :] = param
            errors[i] = error#[sum(abs.(error[1])), sum(abs.(error[2][1]))]
        end
        return collect(params), collect(errors)
    end



    function make_plot1(X, Y, ylabel; 
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf), label=permutedims(L),
        position=:topright, log_scale=:identity, title="")
        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=:identity,
            yaxis=log_scale, framestyle=:box)
        #plot!([], [], label=L"L", alpha=0)

        scatter!(X, Y, label=label, linestyle=:solid,
        legend=position, markershape=markershape, markersize=7,
        markerstrokecolor=colors, framestyle="")
        extra_plot_nr=2
        return p
    end

    extra_plot_nr = 2
    function add_plot(X, Y; label="")
        scatter!(X, Y, label=label, linestyle=:solid, markersize=7,
        markerstrokecolor=colors[extra_plot_nr], framestyle="", markershape=markershape[extra_plot_nr])
        extra_plot_nr += 1
        return p
    end
    labels = permutedims(NR)

    σ_c, x = get_data(L, nr, ts, dist, "most_stressed_fiber",
        "most_stressed_fiber", argmax, ex=[0, 0], rel_x=true, average=false, return_x=true,
        data_path=data_path)
     
 

    LLS_σ_c_N_plot = make_plot1(ts, x[:, :, 1], label=L"σ_c",
        L"k_c/N", log_scale=:identity, title="LLS",
        xlabel=L"t_0", position=:topleft,)

#=     most_stressed_fiber = get_data_kN(L, nr, ts, dist, "average_most_stressed_fiber", return_kN=false)
    E = [σ_to_energy(most_stressed_fiber[1][:, i, 1]) for i in eachindex(ts)]
    k = [argmax(E[1])/length(E[1]) for E in E]
    add_plot(ts, k, label="Max energy") =#


    clusterSize = get_data_kN(L, nr, ts, dist, "average_largest_cluster", return_kN=false)
    localization = [find_localization(clusterSize[1][:, i, 1])/length(clusterSize[1][:, i, 1]) for i in eachindex(ts)]
    add_plot(ts, localization, label="Localization")

    nr_clusters = get_data_kN(L, nr, ts, dist, "average_nr_clusters", return_kN=false)
    localization = [find_localization_nr_clusters(nr_clusters[1][:, i, 1])/length(nr_clusters[1][:, i, 1]) for i in eachindex(ts)]
    add_plot(ts, localization, label="nr_cluster Localization")

    CLS_σ_c_N_plot = make_plot1(ts, x[:, :, 2], label=L"σ_c",
        L"k_c/N", log_scale=:identity, title="CLS", 
        xlabel=L"t_0", position=:topleft,)
    
#=     E = [σ_to_energy(most_stressed_fiber[1][:, i, 2]) for i in eachindex(ts)]
    k = [argmax(E[1])/length(E[1]) for E in E]
    add_plot(ts, k, label="Max energy") =#

    localization = [find_localization(clusterSize[1][:, i, 2])/length(clusterSize[1][:, i, 2]) for i in eachindex(ts)]
    add_plot(ts, localization, label="Localization")

    localization = [find_localization_nr_clusters(nr_clusters[1][:, i, 2])/length(nr_clusters[1][:, i, 2]) for i in eachindex(ts)]
    add_plot(ts, localization, label="nr_cluster Localization")

    return [LLS_σ_c_N_plot,CLS_σ_c_N_plot]
end

L = [128]
α = 2.0
nr = ["LLS", "CLS"]

#ts = vcat((0:20) ./ 50, (5:9) ./ 10)
dist = "ConstantAverageUniform"
#ts = (0:7) ./ 10
data_path = "newData/"
ts = vcat(0.3:0.01:0.5)
ts = vcat(0.05:0.05:0.20, 0.25:0.01:0.5)
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr, dist)
xpsize = 250
ypsize = 250
println("plot")
p = plot(plots..., size=(xpsize * 1.1 * length(plots), ypsize *
                                                       maximum([length(plots) / 2, 1])))

#p2 = plot(plots[3:4]..., size=(psize*length(nr)*1.1,psize
#length(plots)/2/length(nr)), layout = @layout([ A B;]))
display(p)
savefig(p, "plots/Graphs/CriticalStressOverk.pdf")
#savefig(p2, "plots/Graphs/$(dist)_s_over_sigma.pdf")

println("Saved plot!")

