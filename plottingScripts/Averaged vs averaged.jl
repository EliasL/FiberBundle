using Plots
using JLD2
using LaTeXStrings
using Trapz
using LsqFit

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")

funk(X, p) = p[2] .+ p[1] .*X

function get_fit(x, y)
    p0 = [0.22, 0.22]
    return curve_fit(funk, x, y , p0)
end

function lin(x)
    return minimum(x):0.001:maximum(x)
end


function myfit(x, y; fit_interval=1)
    r = round(Int64, fit_interval*length(x))
    println(x[r])
    f = get_fit(x[1:r], y[1:r])
    return funk(lin(x), f.param), f.param
end


function otherPropertiesPlot(L, ts, NR, dist; use_y_lable=true, add_ELS=true)
    
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel, labels; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, series_annotation=[], log=:identity)
        p = scatter(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims, xaxis=log, yaxis=log,
            markershape=[:diamond :rect :star4 :utriangle :dtriangle :circle],
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title, mc=:auto, msc=:auto)
        return p
    end

    function make_plot2(X, Y, ylabel, NR; labels=labels, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, log=:identity)

        p = plot(xlims=xlims, ylims=ylims,  markersize=5, 
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=log, yaxis=log)

        NR = NR=="LLS" ? 1 : 2

        for i in eachindex(X)
            plot!(X[i][:, :, NR], Y[i][:, :, NR], label = labels)
        end
        return p
    end

    function make_plot3(X, Y, ylabel, NR; labels=labels, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf), scale_x=true, scale_y = false,
        position=:topright, log=:identity)

        p = plot(xlims=xlims, ylims=ylims,  markersize=5, 
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=:identity, yaxis=log)

        
        colors = theme_palette(:auto)[1:length(L)]
        markershape=[:diamond :rect :star4 :utriangle :dtriangle :circle]
        series_annotation = [0.3, 0.25]
        series_annotation = text.([t in series_annotation ?  L" $t_0=$"*"$t " : "" for t in ts], pointsize=8, halign=:left, valign=:bottom)

        for i in eachindex(L)
            for nr in NR
                nri = nr=="LLS" ? 1 : 2
                scatter!((scale_x ? X[:, i, nri]./L[i]^2 : X), (scale_y ? Y[:, i, nri]./L[i]^2 : Y[:, i, nri]), label = "$nr "*latexstring("L=$(L[i])"),
                legend=position, markerstrokecolor=colors[i], markershape=markershape[i], series_annotation= (i==1 ? series_annotation : ""))
            end
        end
        return p
    end

    function add_fit!(x, y, fit_interval=0.8)
        y1 = y[:, end, 1]
        f, param = myfit(x, y1, fit_interval=fit_interval)
        param = round.(param, digits=2)
        xx = lin(x)
        s = param[1]>0 ? "+" : "-"
        plot!(xx, f, labels="$(param[2])$(s)$(abs(param[1]))×"*L"t_0", color=:black, linestyle=:dash, alpha=0.5)
    end


    labels = permutedims(NR)
    
    σ_c, x = get_data(L, nr, ts, dist, "most_stressed_fiber", "most_stressed_fiber", argmax, ex=[0,0], average=false, return_x=true, data_path=data_path)
    #σ_c -= [(1-t) / 2 for t=ts, l=L, n = nr]


    σ_c_plot = make_plot3(x, σ_c[:, :, :], log=:identity, 
    L"<σ_c>", permutedims(["$nr" for nr in NR]), title="",
                        xlabel=L"<k_c/N>", position=:topleft, )

    x_c_plot = make_plot3(ts, x, log=:identity, scale_x=false, scale_y=true,
    L"<k_c/N>", permutedims(["$nr" for nr in NR]), title="",
                        xlabel=L"t_0", position=:topleft, )

    σ_cofσ, x = get_data(L, nr, ts, dist, "average_most_stressed_fiber", "most_stressed_fiber", argmax, ex=[0,0], average=true, return_x=true, data_path=data_path)
    #σ_c -= [(1-t) / 2 for t=ts, l=L, n = nr]

    σ_cofσ_plot = make_plot3(x, σ_cofσ[:, :, :], log=:identity, 
    L"σ_c(<σ>)", permutedims(["$nr" for nr in NR]), title="",
                        xlabel=L"k_c/N", position=:topleft, )   

    x_cofσ_plot = make_plot3(ts, x, log=:identity, scale_x=false, scale_y=true,
    L"k_c/N", permutedims(["$nr" for nr in NR]), title="",
                        xlabel=L"t_0", position=:topleft, )
    
    return [σ_c_plot, σ_cofσ_plot, x_c_plot, x_cofσ_plot]
end

L = [16, 32, 64, 128, 256]
α = 2.0
nr = ["LLS"]
ts = vcat(0.05:0.05:0.25, 0.3:0.01:0.5)
dist = "ConstantAverageUniform"
data_path = "newData/"
#ts = (0:7) ./ 10
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr, dist)
xpsize=270
ypsize=330
p = plot(plots..., size=(xpsize*1.1*length(plots)/2,ypsize*length(plots)/2), layout = @layout([ A B ; C D]))
#p2 = plot(plots[3:4]..., size=(psize*length(nr)*1.1,psize*length(plots)/2/length(nr)), layout = @layout([ A B;]))
savefig(p, "plots/Graphs/average or average.pdf")
#savefig(p2, "plots/Graphs/$(dist)_s_over_sigma.pdf")

println("Saved plot!")

