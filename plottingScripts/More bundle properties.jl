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
    upperbounds = [1, 1, -100 , -0.40]
    lowerbounds = [0, 0, -1000, -0.70]
    p0 = [0.1, 0.9, -100, -0.55]
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
        elseif series_position==:right
            sph=:top
        elseif length(series_position)==2
            sph = series_position[2]
            series_position = series_position[1]
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
    div = :max
    largest_cluster_spanning = get_data_spanning(L, nr, ts, "average_largest_cluster", divide=div)
    largest_perimeter_spanning = get_data_spanning(L, nr, ts, "average_largest_perimiter", divide=div)
    most_stressed_fiber_spanning = get_data_spanning(L, nr, ts, "average_most_stressed_fiber", divide=div)
    size_over_σ_LLS = make_plot(largest_cluster_spanning[:, :, 1], log=:log, series_annotation=[[],[],[], [0.14, 0.34]], 
                        L"\tilde{s}_{\mathrm{max}}", permutedims([L"L="*"$l LLS" for l in L]), series_position=:right,
                        x=most_stressed_fiber_spanning[:, :, 1], #= xlabel=L"\tilde{σ}", =# position=:topright, )
    x = most_stressed_fiber_spanning[:, :, 1]
    y = largest_cluster_spanning[:, :, 1]
    f = myfit.(eachcol(x), eachcol(y))
    df = derivative.(eachcol(x), eachcol(y))
    xx = lin.(eachcol(x))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.5,0.35,0.45,0.5)), subplot = 2, xaxis=:log, labels="", xlims=(0.535,0.585),
        xlabel=L"\tilde{σ}", ylabel=L"∂\tilde{s}_{\mathrm{max}} / ∂\tilde{σ}}")

    size_over_σ_CLS = make_plot(largest_cluster_spanning[:, :, 2], log=:log, series_annotation=[[0.7],[],[], [0.34, 0.5]],
                        "",#= L"\tilde{s}_{\mathrm{max}}", =# permutedims([L"L="*"$l CLS" for l in L]), xlims=(-Inf, 1.1),
                        x=most_stressed_fiber_spanning[:, :, 2], #= xlabel=L"\tilde{σ}", =# position=:bottomleft, series_position=:right)
    
    span_over_σ_LLS = make_plot(largest_perimeter_spanning[:, :, 1], log=:log, series_annotation=[[0.0],[],[], [0.12]], 
                        L"\tilde{h}_{\mathrm{max}}", permutedims([L"L="*"$l LLS" for l in L]),
                        x=most_stressed_fiber_spanning[:, :, 1], #= xlabel=L"\tilde{σ}", =# position=:topright, )

    y = largest_perimeter_spanning[:, :, 1]
    f = myfit.(eachcol(x), eachcol(y))
    df = derivative.(eachcol(x), eachcol(y))
    xx = lin.(eachcol(x))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.5,0.35,0.45,0.5)), subplot = 2, xaxis=:log, labels="", xlims=(0.535,0.585),
        xlabel=L"\tilde{σ}", ylabel=L"∂\tilde{h}_{\mathrm{max}} / ∂\tilde{σ}}")

    span_over_σ_CLS = make_plot(largest_perimeter_spanning[:, :, 2], log=:log, series_annotation=[[],[],[],[0.34, 0.6]],
                        "", #= L"\tilde{h}_{\mathrm{max}}",  =#permutedims([L"L="*"$l CLS" for l in L]),
                        x=most_stressed_fiber_spanning[:, :, 2], #= xlabel=L"\tilde{σ}", =# position=:topleft,
                        series_position=[:right, :bottom])
    
    plot!([minimum(most_stressed_fiber_spanning[:, :, 2]),1], [minimum(most_stressed_fiber_spanning[:, :, 2]),1], labels="y=x", color=:black, linestyle=:dash)
    
    x1 = most_stressed_fiber_spanning[:, :, 1]
    x2 = most_stressed_fiber_spanning[:, :, 2]
    y1 = largest_cluster_spanning[:, :, 1] ./ largest_perimeter_spanning[:, :, 1]
    y2 = largest_cluster_spanning[:, :, 2] ./ largest_perimeter_spanning[:, :, 2]
    
    ratio_over_σ_LLS = make_plot(y1, log=:log, series_annotation=[[],[],[],[0.12, 0.30]],
                        L"\tilde{s}_{\mathrm{max}} / \tilde{h}_{\mathrm{max}}", permutedims([L"L="*"$l" for l in L]),
                        x=x1, xlabel=L"\tilde{σ}", position=:topright, title="LLS",
                        series_position=[:right, :bottom])

    f = myfit.(eachcol(x1), eachcol(y1))
    df = derivative.(eachcol(x1), eachcol(y1))
    xx = lin.(eachcol(x1))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.6,0.37,0.36,0.45)), subplot = 2, xaxis=:log, labels="", xlims=(0.535,0.585),
        xlabel=L"\tilde{σ}", ylabel=L"∂\tilde{h}_{\mathrm{max}} / ∂\tilde{σ}}")
    
    ratio_over_σ_CLS = make_plot(y2, log=:log, series_annotation=[[],[],[],[0.34, 0.6]],
                        "", #= L"\tilde{h}_{\mathrm{max}}",  =#permutedims([L"L="*"$l" for l in L]),
                        x=x2, xlabel=L"\tilde{σ}", position=:bottomleft, title="CLS",
                        series_position=[:right, :bottom])
    
    other_plots = [#= size_over_σ_LLS, size_over_σ_CLS, span_over_σ_LLS, span_over_σ_CLS,  =#ratio_over_σ_LLS, ratio_over_σ_CLS]
    return other_plots
end

L = [16, 32, 64, 128]
α = 2.0
nr = ["LLS", "CLS"]
ts = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts2 = (0:9) ./ 10
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr)
p = plot(plots..., size=(620,310), layout = @layout([ A B#= ; C D; E F =#]))
savefig(p, "plots/Graphs/otherBundleProperties.svg")

println("Saved plot!")
