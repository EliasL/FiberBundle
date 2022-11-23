using Plots, Measurements
using JLD2
using LaTeXStrings
using LsqFit


include("../support/ploting_settings.jl")
include("../support/dataManager.jl")

f_lin(x,p) = p[1] .+ x .* p[2]

function get_exp_fit(x,y)

    p0_lin = [0.0, 1.0]
    p0_exp = [1.0]
    fit_lin = curve_fit(f_lin, log.(x), log.(y), p0_lin).param
    f_exp(x,p) = p[1] .* x .^ fit_lin[2]

    fit_exp = curve_fit(f_exp, x, y, p0_exp).param

    return [fit_exp[1], fit_lin[2]]
end

function get_fit(x,y)
    p0_lin = [0.0, 1.0]
    return curve_fit(f_lin, x, y, p0_lin).param
end

function plot_dimension_thing()
    
    L = [8,16,32,64,128,256,512,1024]
    t = [0.0]
    nr = ["SNR", "UNR"]
    N = L.*L
    α = 2.0

    files = [load_file(l, α, t[1], nr[2]) for l in L]
    #seeds = last(files)["nr_seeds_used"]
    #println(seeds)
    s = [f["average_spanning_cluster_size"] for f in files]
    h = [f["average_spanning_cluster_perimiter"] for f in files]

    std_s = [f["std_spanning_cluster_size"] for f in files]
    std_h = [f["std_spanning_cluster_perimiter"] for f in files]

    L = log2.(L)
    s = log2.(s.± std_s)
    h = log2.(h.± std_s)

    fit_s = get_fit(L, Measurements.value.(s))
    fit_h = get_fit(L, Measurements.value.(h))
    slope_s = round(fit_s[2], digits=2)
    slope_h = round(fit_h[2], digits=2)


    s = s .± log2.(std_s)
    h = h .± log2.(std_h)
    m_shape = :cross
    default(color=:black, markersize=3)
    s_plot = scatter(L, s, markershape=m_shape, label=L"Cluster size ($s$)", legend=:topleft, xlabel=L"log$_2(L)$", ylabel=L"log$_2(s)$")
    plot!(L, f_lin(L, fit_s), label=L"Fit ($D_s$)", color=:black)
    h_plot = scatter(L, h, markershape=m_shape, label=L"Perimiter length ($h$)", legend=:topleft, xlabel=L"log$_2(L)$", ylabel=L"log$_2(h)$")
    plot!(L, f_lin(L, fit_h), label=L"Fit ($D_h$)")


    l = @layout [
        A B
    ]
    plot(s_plot, h_plot, size=(500, 300), layout = l,
    plot_title=latexstring("Dimensionality: \$D_s=$slope_s, D_h=$slope_h\$"))

    savefig("plots/Graphs/dimension.pdf")
end

plot_dimension_thing()
println("Saved plot!")