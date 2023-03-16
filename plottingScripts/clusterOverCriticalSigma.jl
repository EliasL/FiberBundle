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
    upperbounds = [ 1,  1, -100 , -0.2]
    lowerbounds = [0.2, 0, -2000, -0.3]
    p0 = [0.2, 0.9, -100, -0.25]
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

function myfunk(x, exp, div)
    return @. x^(exp)/div
end
function otherPropertiesPlot(L, ts, NR, dist; use_y_lable=true, add_ELS=true)
    
    
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
            markershape=[:utriangle :dtriangle :circle :diamond :star4 :rect],
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


    #sigma_C over idssorder
    #smax_over sigma sigma


    labels = permutedims(NR)
    a = [2, 0]
    b = [2, 0]
    c = [0, 0]
    d = [0,0]
    tlabel1=0.3
    tlabel2=0.5
    largest_cluster_spanning = get_data(L, nr, ts, dist, "largest_cluster", "most_stressed_fiber", argmax, ex=a, average=false, data_path=data_path)
    largest_perimeter_spanning = get_data(L, nr, ts, dist, "largest_perimiter", "most_stressed_fiber", argmax, ex=b, average=false, data_path=data_path)
    most_stressed_fiber_spanning = get_data(L, nr, ts, dist, "most_stressed_fiber", "most_stressed_fiber", argmax, ex=c, average=false, data_path=data_path)
    xLLS = most_stressed_fiber_spanning[:, :, 1]
    xCLS = most_stressed_fiber_spanning[:, :, 2]
    xxLLS = lin.(eachcol(xLLS))[end]
    xxCLS = lin.(eachcol(xCLS))[end]
    size_over_σ_LLS = make_plot(largest_cluster_spanning[:, :, 1], log=:log, series_annotation=[[tlabel1, tlabel2],[],[],[], [], []], 
    L"s_{\mathrm{max}}/L^2", permutedims([L"L="*"$l" for l in L]), series_position=:left, title="LLS",
                        x=xLLS, #= xlabel=L"\tilde{σ}", =# position=:bottomleft, )
    
    #plot!(xxLLS[130:end], myfunk(xxLLS[130:end], -8, 2*10^7), labels=L"σ_c^{-8}/(2\times10^7)", color=:black, linestyle=:dash, alpha=0.5)

#=     x = most_stressed_fiber_spanning[:, :, 1]
    y = largest_cluster_spanning[:, :, 1]
    f = myfit.(eachcol(x), eachcol(y))
    df = derivative.(eachcol(x), eachcol(y))
    xx = lin.(eachcol(x))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.5,0.35,0.45,0.5)), subplot = 2, xaxis=:log, labels="", xlims=(xmin, xmax),
        xlabel=L"\tilde{σ}", ylabel=L"∂\tilde{s}_{\mathrm{max}} / ∂\tilde{σ}}") =#

    size_over_σ_CLS = make_plot(largest_cluster_spanning[:, :, 2], log=:log, series_annotation=[[tlabel1, tlabel2],[],[],[], [], []],
    L"s_{\mathrm{max}}", permutedims([L"L="*"$l" for l in L]), title="CLS",
                        x=xCLS, #= xlabel=L"\tilde{σ}", =# position=:bottomleft, series_position=:left)
    
    #plot!(xxCLS, myfunk(xxCLS, -5.5, 50), labels=L"σ_c^{-5.5}/50", color=:black, linestyle=:dash, alpha=0.5)


    span_over_σ_LLS = make_plot(largest_perimeter_spanning[:, :, 1], log=:log, series_annotation=[[tlabel1, tlabel2],[],[],[], [], []], 
    L"h_{\mathrm{max}}/L^2", permutedims([L"L="*"$l" for l in L]),
                        x=xLLS, #= xlabel=L"\tilde{σ}", =# position=:bottomleft, )
                            
    #plot!(xxLLS[100:end], myfunk(xxLLS[100:end], -7, 5*10^6), labels=L"σ_c^{-6}/(4\times10^5)", color=:black, linestyle=:dash, alpha=0.5)

#=     y = largest_perimeter_spanning[:, :, 1]
    f = myfit.(eachcol(x), eachcol(y))
    df = derivative.(eachcol(x), eachcol(y))
    xx = lin.(eachcol(x))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.5,0.35,0.45,0.5)), subplot = 2, xaxis=:log, labels="", xlims=(xmin, xmax),
        xlabel=L"\tilde{σ}", ylabel=L"∂\tilde{h}_{\mathrm{max}} / ∂\tilde{σ}}") =#
       

    span_over_σ_CLS = make_plot(largest_perimeter_spanning[:, :, 2], log=:log, series_annotation=[[tlabel1, tlabel2],[],[],[], [], []],
    L"h_{\mathrm{max}}", permutedims([L"L="*"$l" for l in L]),
                        x=xCLS, #= xlabel=L"\tilde{σ}", =# position=:bottomleft,
                        series_position=:left)
    
    #plot!(xxCLS, myfunk(xxCLS, -4, 6), labels=L"σ_c^{-4}/6", color=:black, linestyle=:dash, alpha=0.5)
    #plot!([minimum(xCLS),maximum(xCLS)], [minimum(xCLS),maximum(xCLS)], labels="y=x", color=:black, linestyle=:dash)

    largest_cluster_spanning = get_data(L, nr, ts, dist, "largest_cluster", "most_stressed_fiber", argmax, ex=[0,0], average=false, data_path=data_path)
    largest_perimeter_spanning = get_data(L, nr, ts, dist, "largest_perimiter", "most_stressed_fiber", argmax, ex=[0,0], average=false, data_path=data_path)

    y1 = (largest_cluster_spanning[:, :, 1] ./ largest_perimeter_spanning[:, :, 1]) ./ permutedims(L.^d[1])
    y2 = (largest_cluster_spanning[:, :, 2] ./ largest_perimeter_spanning[:, :, 2]) ./ permutedims(L.^d[2])
    
    ratio_over_σ_LLS = make_plot(y1, log=:log, series_annotation=[[tlabel1, tlabel2],[],[],[], [], []],
                        L" s_{\mathrm{max}}/h_{\mathrm{max}}", permutedims([L"L="*"$l" for l in L]),
                        x=xLLS, xlabel=L"σ_c", position=:bottomleft, title="LLS",
                        series_position=:left)
    
    #plot!(xxLLS, myfunk(xxLLS, -3, 30), labels=L"σ_c^{-3}/30", color=:black, linestyle=:dash, alpha=0.5)
    #plot!(xxLLS, myfunk(xxLLS, -1, 3), labels=L"σ_c^{-1}/3", color=:black, linestyle=:dot, alpha=0.5)
#= 
    f = myfit.(eachcol(x1), eachcol(y1))
    df = derivative.(eachcol(x1), eachcol(y1))
    xx = lin.(eachcol(x1))
    plot!(xx, f, labels="", color=permutedims(theme_palette(:auto)[1:length(L)]))
    plot!(xx, df, inset = (1, bbox(0.6,0.5,0.36,0.35)), subplot = 2, xaxis=:log, labels="", xlims=(xmin, xmax),
        xlabel=""#= L"\tilde{σ}" =#, ylabel=L"∂\tilde{h}_{\mathrm{max}} / ∂\tilde{σ}}") =#
    
    ratio_over_σ_CLS = make_plot(y2, log=:log, series_annotation=[[tlabel1, tlabel2],[],[],[], [], []],
    L" s_{\mathrm{max}}/h_{\mathrm{max}}", #= L"\tilde{h}_{\mathrm{max}}",  =#permutedims([L"L="*"$l" for l in L]),
                        x=xCLS, xlabel=L"σ_c", position=:bottomleft, title="CLS",
                        series_position=:left)

    #plot!(xxCLS, myfunk(xxCLS, -1.75, 10), labels=L"σ_c^{-1.75}/10", color=:black, linestyle=:dash, alpha=0.5)
    
    other_plots = [size_over_σ_LLS, size_over_σ_CLS, span_over_σ_LLS, span_over_σ_CLS, ratio_over_σ_LLS, ratio_over_σ_CLS]
    return other_plots
end

L = [32, 64, 128, 256]#[16, 32, 64, 128, 256, 512]
α = 2.0
nr = ["LLS", "CLS"]
dist = "ConstantAverageUniform"
data_path="newData/"
#ts = vcat((0:20) ./ 50, (5:7) ./ 10)
ts = vcat(0.15:0.05:0.25, 0.3:0.01:0.5)
#ts2 = (0:9) ./ 10
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = otherPropertiesPlot(L, ts, nr, dist)
psize=300
p = plot(plots..., size=(psize*length(nr)*1.2,psize*length(plots)/length(nr)), layout = @layout([ A B; C D; E F]))
savefig(p, "plots/Graphs/$(dist)_sh_over_sigma_c.pdf")

println("Saved plot!")
