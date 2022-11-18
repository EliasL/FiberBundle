using Plots
using JLD2
using LaTeXStrings
using LsqFit


include("../support/ploting_settings.jl")
include("../support/dataManager.jl")


f(x,p) = p[1].*x.^p[2]
function get_fit(x,y)
    p0 = [5, 1.7]
    fit = curve_fit(f, x, y, p0)
    return fit.param  
end

function plot_line()
    
end

function plot_dimension_thing()
    
    L = [8,16,32,64,128,256]
    t = [0.0]
    nr = ["SNR", "UNR"]
    N = L.*L
    α = 2.0

    files = [load_file(l, α, t[1], nr[2]) for l in L]
    seeds = last(files)["nr_seeds_used"]
    println(seeds)
    s = [f["average_spanning_cluster_size"] for f in files]
    h = [f["average_spanning_cluster_perimiter"] for f in files]
    
    fit_s = get_fit(L, s)
    fit_h = get_fit(L, h)
    slope_s = round(fit_s[2], digits=2)
    slope_h = round(fit_h[2], digits=2)

    p = plot(L, s, xaxis=:log, yaxis=:log, seriestype = :scatter, label="s", legend=:topleft, xlabel=L"L", ylabel=L"s",
             plot_title=latexstring("Dimensions $(nr[2]) \$t_0 = $(t[1]), D_s=$slope_s, D_h=$slope_h\$"), plot_titlevspan=0.1)
    plot!(L, h, seriestype = :scatter, label="h")
    xticks!(p, L, string.(L))
    plot!(L, f(L, fit_s), label="s fit")
    plot!(L, f(L, fit_h), label="h fit")

    savefig("plots/Graphs/dimension.pdf")
end

plot_dimension_thing()
println("Saved plot!")