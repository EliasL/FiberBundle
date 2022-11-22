using Plots, Measurements
using JLD2
using LaTeXStrings
using LsqFit


include("../support/ploting_settings.jl")
include("../support/dataManager.jl")


function get_fit(x,y)

    f_lin(x,p) = p[1] .+ x .* p[2]
    p0_lin = [0.0, 1.0]
    p0_exp = [1.0]
    fit_lin = curve_fit(f_lin, log.(x), log.(y), p0_lin).param
    f_exp(x,p) = p[1] .* x .^ fit_lin[2]

    fit_exp = curve_fit(f_exp, x, y, p0_exp).param

    return [fit_exp[1], fit_lin[2]]
end

function plot_dimension_thing()
    
    L = [8,16,32,64,128,256,512,1024]
    t = [0.0]
    nr = ["SNR", "UNR"]
    N = L.*L
    α = 2.0

    files = [load_file(l, α, t[1], nr[2]) for l in L]
    seeds = last(files)["nr_seeds_used"]
    println(seeds)
    s = [f["average_spanning_cluster_size"] for f in files]
    h = [f["average_spanning_cluster_perimiter"] for f in files]

    std_s = [f["std_spanning_cluster_size"] for f in files]
    std_h = [f["std_spanning_cluster_perimiter"] for f in files]

    fit_s = get_fit(L, s)
    fit_h = get_fit(L, h)
    println(fit_s)
    slope_s = round(fit_s[2], digits=2)
    slope_h = round(fit_h[2], digits=2)


    s = s .± std_s

    p = plot(L, s, xaxis=:log, yaxis=:log, seriestype = :scatter, label="s", legend=:topleft, xlabel=L"L", ylabel=L"s",
             plot_title=latexstring("Dimensions $(nr[2]) \$t_0 = $(t[1]), D_s=$slope_s"),#, D_h=$slope_h\$"),
             plot_titlevspan=0.1)
    #plot!(L, h, seriestype = :scatter, label="h", ribbon=[(std_s,std_s)])
    xticks!(p, L, string.(L))
    plot!(L, fit_s[1] .*L.^slope_s, label="s fit")
    #plot!(L, f(L, fit_h), label="h fit")
    savefig("plots/Graphs/dimension.pdf")
end

plot_dimension_thing()
println("Saved plot!")