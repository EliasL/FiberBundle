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


function myfit(x, y; fit_interval=1)
    r = round(Int64, fit_interval * length(x))
    f = get_fit(x[1:r], y[1:r])
    return funk(lin(x), f.param), f.param, margin_error(f,)
end

function otherPropertiesPlot(L, ts, NR, dist; use_y_lable=true)
    yLabel(string) = use_y_lable ? string : ""

    colors = permutedims(theme_palette(:auto)[1:16])
    markershape = [:circle :pentagon :star4 :diamond :rect :ltriangle :rtriangle :star6]

    function add_fits(x, y, plot_index)
        params = zeros(Float64, (size(y)[1],2))
        errors = zeros(Float64, (size(y)[1],2))
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
            errors[i, :] = error#[sum(abs.(error[1])), sum(abs.(error[2][1]))]
        end
        return collect(params), collect(errors)
    end



    function make_plot1(X, Y, ylabel; 
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topright, log_scale=:identity)

        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=" ", xaxis=:identity,
            yaxis=log_scale)
        plot!([], [], label=L"t_0", alpha=0)

        plot_ts = [0.1, 0.30, 0.35, 0.4, 0.5]
        plot_index = [i for i in eachindex(ts) if ts[i] in plot_ts] 
        labels = permutedims([latexstring("$t") for t in plot_ts])
            
        plot_Y = permutedims(Y[[i for i in eachindex(ts) if ts[i] in plot_ts], :, 2])
        scatter!(X, plot_Y, label=labels, linestyle=:solid,
        legend=position, markershape=markershape, markersize=7,
        markerstrokecolor=colors)
        
        params, errors = add_fits(X, Y[:, :, 2], plot_index)
        #println((errors[:,1]))
        scatter!(ts, [params[:, 1] errors[:,2] ], label=["Slope" "Error"],
            markershape=markershape, markerstrokecolor=colors, markersize=7,
            legend=:topleft, xlabel=L"t_0", ylabel="",
            inset = (1, bbox(0.40, 0.42, 0.58, 0.4)), subplot=2)
        return p
    end

    labels = permutedims(NR)

    σ_c, x = get_data(L, nr, ts, dist, "most_stressed_fiber",
        "most_stressed_fiber", argmax, ex=[0, 0], average=false, return_x=true,
        data_path=data_path)

    N = L .* L
    lnN =log.(log.(N))
    σ_c_N_plot = make_plot1(lnN, 1 ./ σ_c[:, :, :],
        L"1/ \langle σ_c \rangle", log_scale=:identity,
        xlabel=L"\ln(\ln(N))", position=:bottomleft,)


    return [σ_c_N_plot]
end

L = [16, 32, 64, 128, 256, 512]
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
xpsize = 330
ypsize = 330
p = plot(plots..., size=(xpsize * 1.1 * length(plots), ypsize *
                                                       maximum([length(plots) / 2, 1])))

#p2 = plot(plots[3:4]..., size=(psize*length(nr)*1.1,psize
#length(plots)/2/length(nr)), layout = @layout([ A B;]))
display(p)
savefig(p, "plots/Graphs/CriticalStressOverN.pdf")
#savefig(p2, "plots/Graphs/$(dist)_s_over_sigma.pdf")

println("Saved plot!")

