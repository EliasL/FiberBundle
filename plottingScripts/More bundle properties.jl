using Plots
using JLD2
using LaTeXStrings
using Trapz
using LsqFit

include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("../support/dataIterator.jl")

@. sigmoid(x,p) = p[1] + p[2] / (1 +ℯ^(p[3]*(x+p[4])))
@. dSogmoid(x,p) = -p[2]*p[3]*ℯ^(p[3]*(x +p[4])) / (ℯ^(p[3]*(x +p[4])) + 1)^2

function get_fit(x, y)
    upperbounds = [1, 1, -100 , -0.50]
    lowerbounds = [0, 0, -1000, -0.60]
    p0 = [0.1, 0.9, -200, -0.55]
    return curve_fit(sigmoid, x, y , p0, lower=lowerbounds, upper=upperbounds)
end

function lin(x)
    return minimum(x):0.001:maximum(x)
end

function derivative(x, y)
    f = get_fit(x, y)
    df(x) = dSogmoid(x, f.param)
    return df.(lin(x)) 
end

function myfit(x, y)
    f = get_fit(x, y)
    return sigmoid(lin(x), f.param)
end

function otherPropertiesPlot(L, ts, NR; use_y_lable=true, add_ELS=true)
    
    
    yLabel(string) = use_y_lable ? string : ""
    function make_plot(y, ylabel, labels; x=k_N, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, series_annotation=[], series_position=:left, log=:identity, anotate_both=true)
        if series_position==:left
            sph=:bottom
        else
            sph=:bottom
        end
        if !anotate_both
            series_annotation = [["" for _ in ts] text.([t in series_annotation ?  L" t_0 = "*"$t" : "" for t in ts], halign=series_position, valign=sph)]
        else
            series_annotation = permutedims([text.([t in sa ?  L" $t_0=$"*"$t " : "" for t in ts], halign=series_position, valign=sph, pointsize=8) for sa in series_annotation])
        end
        p = scatter(x, y, label = labels, legend=position, xlims=xlims, ylims=ylims, xaxis=log, yaxis=log,
            markershape=[:utriangle :diamond :circle],
            xlabel=xlabel, ylabel=yLabel(ylabel), title=title, mc=:auto, msc=:auto, 
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
    largest_perimeter_spanning = get_data_spanning(L, nr, ts, "average_largest_perimiter", divide=:max)
    most_stressed_fiber_spanning = get_data_spanning(L, nr, ts, "average_most_stressed_fiber", divide=:max)
    size_over_σ_LLS = make_plot(largest_cluster_spanning[:, :, 1], log=:log, series_annotation=[[0.00],[],[]], 
                        L"\tilde{s}_{\mathrm{max}}", permutedims([L"L="*"$l LLS" for l in L]),
                        x=most_stressed_fiber_spanning[:, :, 1], #= xlabel=L"\tilde{σ}", =# position=:topright, )
    x = most_stressed_fiber_spanning[:, :, 1]
    y = largest_cluster_spanning[:, :, 1]
    f = myfit.(eachcol(x), eachcol(y))
    df = derivative.(eachcol(x), eachcol(y))
    xx = lin.(eachcol(x))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.5,0.35,0.45,0.5)), subplot = 2, xaxis=:log, labels="", xlims=(0.535,0.585),
        xlabel=L"\tilde{σ}", ylabel=L"∂\tilde{s}_{\mathrm{max}} / ∂\tilde{σ}}")

    size_over_σ_CLS = make_plot(largest_cluster_spanning[:, :, 2], log=:log, series_annotation=[[0.7],[],[]],
                        "",#= L"\tilde{s}_{\mathrm{max}}", =# permutedims([L"L="*"$l CLS" for l in L]), xlims=(-Inf, 1.1),
                        x=most_stressed_fiber_spanning[:, :, 2], #= xlabel=L"\tilde{σ}", =# position=:bottomleft, series_position=:left)
    
    span_over_σ_LLS = make_plot(largest_perimeter_spanning[:, :, 1], log=:log, series_annotation=[[0.0],[],[]], 
                        L"\tilde{h}_{\mathrm{max}}", permutedims([L"L="*"$l LLS" for l in L]),
                        x=most_stressed_fiber_spanning[:, :, 1], xlabel=L"\tilde{σ}", position=:topright, )

    y = largest_perimeter_spanning[:, :, 1]
    f = myfit.(eachcol(x), eachcol(y))
    df = derivative.(eachcol(x), eachcol(y))
    xx = lin.(eachcol(x))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.5,0.35,0.45,0.5)), subplot = 2, xaxis=:log, labels="", xlims=(0.535,0.585),
        xlabel=L"\tilde{σ}", ylabel=L"∂\tilde{h}_{\mathrm{max}} / ∂\tilde{σ}}")

    span_over_σ_CLS = make_plot(largest_perimeter_spanning[:, :, 2], log=:identity, series_annotation=[[],[],[0.04]],
                        "", #= L"\tilde{h}_{\mathrm{max}}",  =#permutedims([L"L="*"$l CLS" for l in L]), xlims=(0,1), ylims=(0,1),
                        x=most_stressed_fiber_spanning[:, :, 2], xlabel=L"\tilde{σ}", position=:topleft, series_position=:right)
    
    plot!([0,1], [0,1], labels="y=x", color=:black, linestyle=:dash)
    other_plots = [size_over_σ_LLS, size_over_σ_CLS, span_over_σ_LLS, span_over_σ_CLS]
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
p = plot(plots..., size=(620,620))
savefig(p, "plots/Graphs/otherBundleProperties.svg")

println("Saved plot!")
