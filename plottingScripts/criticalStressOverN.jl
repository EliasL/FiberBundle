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
    return LinRange(minimum(x),maximum(x),10)
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
    markershape = [:circle :pentagon :star4 :diamond :rect :utriangle :dtriangle :star6]

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
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:topleft, log_scale=:identity)

        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=" ", xaxis=:identity,
            yaxis=log_scale, framestyle=:box)
        plot!([], [], label=L"t_0", alpha=0)

        plot_ts = [0.1, 0.25, 0.3, 0.4, 0.5]
        plot_index = [i for i in eachindex(ts) if ts[i] in plot_ts] 
        labels = permutedims([latexstring("$t") for t in plot_ts])
            
        plot_Y = permutedims(Y[[i for i in eachindex(ts) if ts[i] in plot_ts], :, 2])
        scatter!(X, plot_Y, label=labels, linestyle=:solid,
        legend=position, markershape=markershape, markersize=7,
        markerstrokecolor=colors, framestyle="", ylims=(-Inf, 0.41))
        
        params, errors = add_fits(X, Y[:, :, 2], plot_index)
        #println((errors[:,1]))
        scatter!(ts, [params[:, 1] ], label="Slope", framestyle="", 
            markershape=markershape, markerstrokecolor=colors, markersize=7,
            legend=:topleft, xlabel=L"t_0", ylabel="", ylims=(-Inf, 0.35),
            inset = (1, bbox(0.37, 0.1, 0.50, 0.4)), subplot=2)
        scatter!( twinx(p[2]), ts, errors, label="Dev.", legend=:topright,
            framestyle="", ylims=(0, 0.003),
            markershape=markershape[3], markerstrokecolor=colors[2], markersize=7,
        )
        
        
        return p
    end

    function make_plot2(X, Y, ylabel; 
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:bottomleft, log_scale=:identity)

        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=" ", xaxis=:identity,
            yaxis=log_scale, framestyle=:box)
        plot!([], [], label=L"t_0", alpha=0)

        plot_ts = [0.1, 0.2, 0.25, 0.26, 0.5]
        plot_index = [i for i in eachindex(ts) if ts[i] in plot_ts] 
        labels = permutedims([latexstring("$t") for t in plot_ts])
            
        plot_Y = permutedims(Y[[i for i in eachindex(ts) if ts[i] in plot_ts], :, 2])
        scatter!(X, plot_Y, label=labels, linestyle=:solid,
        legend=position, markershape=markershape, markersize=7,
        markerstrokecolor=colors, framestyle="")
        
        params, errors = add_fits(X, Y[:, :, 2], plot_index)
        #println((errors[:,1]))
        scatter!(ts, [params[:, 1] ], label="Slope", framestyle="", 
            markershape=markershape, markerstrokecolor=colors, markersize=7,
            legend=:topleft, xlabel=L"t_0", ylabel="", ylims=(-Inf, 1.3),
            inset = (1, bbox(0.368, 0.32, 0.50, 0.4)), subplot=2)
        scatter!( twinx(p[2]), ts, errors, label="Dev.", legend=:topright,
            framestyle="", ylims=(0, 0.010),
            markershape=markershape[3], markerstrokecolor=colors[2], markersize=7,
        )
        
        
        return p
    end

    function make_plot3(X, Y, ts, ylabel; 
        ylims=(-Inf, Inf), xlabel="", xlims=(-Inf, Inf),
        position=:outerright, log_scale=:identity)

        p = plot(xlims=xlims, ylims=ylims, markersize=5,
            xlabel=xlabel, ylabel=yLabel(ylabel), title=" ", xaxis=:identity,
            yaxis=log_scale, framestyle=:box)
        plot!([], [], label=L"t_0", alpha=0)

        plot_ts = ts
        plot_index = [i for i in eachindex(ts) if ts[i] in plot_ts] 
        labels = permutedims([latexstring("$t") for t in plot_ts])
            
        plot_Y = permutedims(Y[[i for i in eachindex(ts) if ts[i] in plot_ts], :, 2])
        scatter!(X, plot_Y, label=labels, linestyle=:solid,
        legend=position, markershape=markershape, markersize=7,
        markerstrokecolor=colors, framestyle="", ylims=(-0.1, Inf))
        
        params, errors = add_fits(X, Y[:, :, 2], plot_index)
        #println((errors[:,1]))
        scatter!(ts, [params[:, 1] ], label="Slope", framestyle="", 
            markershape=markershape, markerstrokecolor=colors, markersize=7,
            legend=:topright, xlabel=L"t_0", ylabel="", ylims=(-Inf, 0.33),
            inset = (1, bbox(0.5, 0.50, 0.45, 0.35)), subplot=2)
        return p
    end

    labels = permutedims(NR)

    σ_c, k_c = get_data(L, nr, ts, dist, "most_stressed_fiber",
        "most_stressed_fiber", argmax, ex=[0, 0], average=false, return_x=true,
        data_path=data_path, rel_x=true)

    N = L .* L
    lnN = 1 ./ log.(log.(N))
    σ_c_N_plot = make_plot1(lnN, σ_c[:, :, :],
        L"\langle σ_c \rangle", log_scale=:identity,
        xlabel=L"1/\ln(\ln(N))")
    

    k_c_N_plot = make_plot2(lnN, k_c[:, :, :],
        L"\langle k_c \rangle", log_scale=:identity,
        xlabel=L"1/\ln(\ln(N))")

    tts = [0.1, 0.15, 0.20, 0.25, 0.26, 0.27, 0.28, 0.29, 0.31, 0.32, 0.33, 0.34, 0.35, 0.40, 0.45, 0.5]
    σ_c, k_c = get_data(L, nr, tts, dist, "most_stressed_fiber",
    "most_stressed_fiber", argmax, ex=[0, 0], average=false, return_x=true,
    data_path=data_path, rel_x=true)
#=     r_k_c =(reverse(k_c, dims=1))
    strange_p = make_plot3(lnN, σ_c .- r_k_c, tts,
        L"\langle \sigma_c \rangle_{t_0} - \langle k_c \rangle_{0.6-t_0} ", log_scale=:identity,
        xlabel=L"1/\ln(\ln(N))") =#
    r_k_c =(reverse(k_c, dims=1))
    println(σ_c[:, 1, 2])
    println(r_k_c[:, 1, 2])
    strange_p = scatter(tts, σ_c[:, 6, 2] .- r_k_c[:, 6, 2])
    return [σ_c_N_plot, k_c_N_plot, strange_p]
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
σ_plot, k_plot, strange_plot = otherPropertiesPlot(L, ts, nr, dist)
xpsize = 360
ypsize = 330
p = plot(σ_plot, size=(xpsize * 1.1 * length(plots), ypsize *
                                                       maximum([length(plots) / 2, 1])))
p2 = plot(k_plot, size=(xpsize * 1.1 * length(plots), ypsize *
                                                       maximum([length(plots) / 2, 1])))

p3 = plot(strange_plot, size=(400, 350))
#p2 = plot(plots[3:4]..., size=(psize*length(nr)*1.1,psize
#length(plots)/2/length(nr)), layout = @layout([ A B;]))
#display(p2)
savefig(p, "plots/Graphs/CriticalStressOverN.pdf")
savefig(p2, "plots/Graphs/CriticalStressKOverN.pdf")
savefig(p3, "plots/Graphs/StrangeCriticalStressAndKOverN.pdf")
#savefig(p2, "plots/Graphs/$(dist)_s_over_sigma.pdf")

println("Saved plot!")

