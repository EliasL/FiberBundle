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

    function make_plot3(X, Y, ylabel, NR; labels=labels, title="", ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, log=:identity)

        p = plot(xlims=xlims, ylims=ylims,  markersize=5, 
        xlabel=xlabel, ylabel=yLabel(ylabel), title=title, xaxis=:identity, yaxis=:identity)

        
        colors = theme_palette(:auto)[1:length(NR)]
        markershape=[:diamond :rect :star4 :utriangle :dtriangle :circle]

        for i in eachindex(L)
            for nr in NR
                nri = nr=="LLS" ? 1 : 2
                scatter!(X, Y[:, i, nri], label = "$nr "*latexstring("L=$(L[i])"),
                legend=position, markerstrokecolor=colors[nri], markershape=markershape[i])
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
    
        #= y2 = y[:, end, 2]
        f, param = myfit(x, y2, fit_interval=fit_interval)
        param = round.(param, digits=2)
        xx = lin(x)
        plot!(xx, f, labels="$(param[2])+$(param[1])×"*L"t_0", color=:black, linestyle=:dot, alpha=0.5)
         =#
    end

    #sigma_C over idssorder
    #smax_over sigma sigma


    labels = permutedims(NR)
    

    σ_c = get_data(L, nr, ts, dist, "most_stressed_fiber", "most_stressed_fiber", argmax, ex=[0,0], average=false, data_path=data_path)
    #σ_c -= [(1-t) / 2 for t=ts, l=L, n = nr]

    ts = vcat([0.0], ts)
    new_σ_c = zeros(length(ts), length(L), length(nr))
    for l=eachindex(L), n=eachindex(nr)
        new_σ_c[:, l, n] = vcat([0.5], σ_c[:, l, n]) 
    end

    σ_c_plot = make_plot3(ts, new_σ_c, log=:log, 
    L"\langle σ_c \rangle", permutedims(["$nr" for nr in NR]), title="",
                        xlabel=xlabel, position=:topright, )

    #add_fit!(ts, σ_c)
#= 
    σ_cofσ = get_data(L, nr, ts, dist, "average_most_stressed_fiber", "most_stressed_fiber", argmax, ex=[0,0], average=true, data_path=data_path)
    #σ_c -= [(1-t) / 2 for t=ts, l=L, n = nr]
    σ_cofσ_plot = make_plot3(ts, σ_cofσ[:, :, :], log=:log, 
    L"σ_c(\langle σ \rangle)", permutedims(["$nr" for nr in NR]), title="",
                        xlabel=L"t_0", position=:topright, ) =#
    #add_fit!(ts, σ_c)
    max_y = 0.5#round(maximum(σ_c[:, :, 1]), digits=2)
    slope = 1
    #a = 22
    #b = 0.54
    #plot!(ts[1:6], max_y .- slope.*ts[1:6], labels="$max_y"*L" -" * L"t_0", color=:black, linestyle=:dash, alpha=0.5)

    
    
    return [σ_c_plot]
end

α = 2.0
nr = ["LLS", "CLS"]

dist = "ConstantAverageUniform"
dist = "Weibull"
if dist=="Weibull"
    L=[128]
    xlabel=L"t_w"
    ts = vcat([0.5], 1:0.1:1.5, 2:0.5:5)
else
    L = [64, 128, 256, 512]
    pos = :topleft
    xlabel=L"t_0"
    ts = vcat(0.05:0.05:0.20, 0.25:0.01:0.5)
end
data_path="newData/"
#ts = (0:7) ./ 10
#ts2 = vcat((0:20) ./ 50, (5:9) ./ 10)
#ts = [0.1,0.2]
plots = [otherPropertiesPlot(L, ts, nr, dist)[1]]
xpsize=300
ypsize=300
p = plot(plots..., size=(xpsize*1.1*length(plots),ypsize), layout = @layout([ A B ;]))
#p2 = plot(plots[3:4]..., size=(psize*length(nr)*1.1,psize*length(plots)/2/length(nr)), layout = @layout([ A B;]))
savefig(p, "plots/Graphs/$(dist)_average_sigma_C_over_t0.pdf")
#savefig(p2, "plots/Graphs/$(dist)_s_over_sigma.pdf")

println("Saved plot!")

