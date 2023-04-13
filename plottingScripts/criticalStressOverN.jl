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
    return minimum(x):0.001:maximum(x)
end


function myfit(x, y; fit_interval=1)
    r = round(Int64, fit_interval * length(x))
    println(x[r])
    f = get_fit(x[1:r], y[1:r])
    return funk(lin(x), f.param), f.param
end


function otherPropertiesPlot(L, ts, NR, dist; use_y_lable=true, add_ELS=true)


    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel, labels; x=k_N, title="", ylims=(-Inf, Inf),
        xlabel="", xlims=(-Inf, Inf),
        position=:topright, series_annotation=[], log=:identity)
        p = scatter(x, y, label=labels, legend=position, xlims=xlims,
            ylims=ylims, xaxis=log, yaxis=log,
            markershape=[:diamond :rect :star4 :utriangle :dtriangle :circle],
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title, mc=:auto,
            msc=:auto)
        return p
    end

    function make_plot2(X, Y, ylabel, NR; labels=labels, title="",
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, log=:identity)

        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=log,
            yaxis=log)

        NR = NR == "LLS" ? 1 : 2

        for i in eachindex(X)
            plot!(X[i][:, :, NR], Y[i][:, :, NR], label=labels)
        end
        return p
    end

    function make_plot3(X, Y, ylabel, NR; labels=labels, title="",
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, log_scale=:identity)

        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=:identity,
            yaxis=log_scale)


        colors = permutedims(vcat(theme_palette(:auto), theme_palette(:auto)))
        markershape = [:diamond :rect :star4 :utriangle :dtriangle :circle]
        series_annotation = [0.3, 0.5]
        series_annotation = text.([t in series_annotation ? L" $t_0=$"*"$t " : "" for t in ts], pointsize=8,
            halign=:left, valign=:bottom)
        labels = permutedims([latexstring("t_0=$(ts[i])") for i in eachindex(ts)])
        plot!(X, permutedims(Y), label=labels,
        markerstrokecolor=permutedims(theme_palette(:auto)[1:length(ts)]),
        legend=position, markershape=markershape)
        plot!(X, X.+3.2, label=L"\ln(\ln(N))+3.2", c=:black)

        return p
    end

    function make_plot4(X, Y, ylabel; 
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, log_scale=:identity)

        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=" ", xaxis=:identity,
            yaxis=log_scale)


        colors =theme_palette(:auto)
        markershape = [:diamond :rect :star4 :utriangle :dtriangle :circle]

        series_annotation = [0.3, 0.5]
        series_annotation = text.([t in series_annotation ? L" $t_0=$"*"$t " : "" for t in ts], pointsize=8,
            halign=:left, valign=:bottom)
        labels = permutedims([latexstring("t_0=$(ts[i])") for i in eachindex(ts)])
         
        plot!([],[[],[],[],[],[]], label=labels,
        c=:black, markerstrokecolor=:black,
        legend=position, markershape=markershape)
         # Add coloured text in the middle of the plot; each subsequent group name
        # is shifted by an amount that is proportional to maximum(length.(group_names))
        
        plot!(X, permutedims(Y[:, :, 2]), label="", linestyle=:solid,
        c=colors[2], markerstrokecolor=colors[2],
        legend=position, markershape=markershape)       
        plot!(X, permutedims(Y[:, :, 1]), label="",
        c=colors[1], markerstrokecolor=colors[1], linestyle=:dash,
        legend=position, markershape=markershape)
        plot!(X, 1.25*X.+2.1, label=L"\frac{5}{4}\ln(\ln(N))+2.1", c=:black)
        plot!([], [], label=" ", alpha=0)
        tc = [colors[1], :black, colors[2]]
        for (i,s) in enumerate(["LLS", "and", "CLS"])
            annotate!(1.85+ 0.135*i , 5.7, text(s, tc[i], 14))
        end
        return p
    end

    labels = permutedims(NR)

    σ_c, x = get_data(L, nr, ts, dist, "most_stressed_fiber",
        "most_stressed_fiber", argmax, ex=[0, 0], average=false, return_x=true,
        data_path=data_path)
    #σ_c -= [(1-t) / 2 for t=ts, l=L, n = nr]

    N = L .* L
    lnN = log.(log.(N))
#=     LLS_σ_c_N_plot = make_plot3(lnN, 1 ./ σ_c[:, :, 1],
        L"1/\langle σ_c \rangle", "LLS",
        labels=permutedims(["$nr" for nr in NR]), title="LLS", log_scale=:identity,
        xlabel=L"\ln(N^2)/\ln(\ln(N^2))", position=:bottomright,)

    CLS_σ_c_N_plot = make_plot3(lnN, 1 ./ σ_c[:, :, 2],
        L"1/ \langle σ_c \rangle", "CLS",
        labels=permutedims(["$nr" for nr in NR]), title="CLS", log_scale=:identity,
        xlabel=L"\ln(N^2)/\ln(\ln(N^2))", position=:bottomright,)
 =#    
    σ_c_N_plot = make_plot4(lnN, 1 ./ σ_c[:, :, :],
        L"1/ \langle σ_c \rangle", log_scale=:identity,
        xlabel=L"\ln(\ln(N))", position=:bottomright,)


    return [σ_c_N_plot]#[LLS_σ_c_N_plot, CLS_σ_c_N_plot]
end

L = [16, 32, 64, 128, 256, 512]
α = 2.0
nr = ["LLS", "CLS"]

#ts = vcat((0:20) ./ 50, (5:9) ./ 10)
dist = "ConstantAverageUniform"
#ts = (0:7) ./ 10
data_path = "newData/"
ts = reverse!([0.05, 0.1, 0.3, 0.35, 0.5])
#ts = vcat(0.05:0.05:0.25, 0.3:0.01:0.5)
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr, dist)
xpsize = 270
ypsize = 330
p = plot(plots..., size=(xpsize * 1.1 * length(plots), ypsize *
                                                       maximum([length(plots) / 2, 1])))
#p2 = plot(plots[3:4]..., size=(psize*length(nr)*1.1,psize
#length(plots)/2/length(nr)), layout = @layout([ A B;]))
savefig(p, "plots/Graphs/CriticalStressOverN.pdf")
#savefig(p2, "plots/Graphs/$(dist)_s_over_sigma.pdf")

println("Saved plot!")
p

