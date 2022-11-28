using Plots, Measurements
using JLD2
using LaTeXStrings
using LsqFit


include("../support/ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/inertia.jl")

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
    
    L = [8,16,32,64,128,256,512]
    t = (0:9)./10
    NR = ["SNR", "UNR"]
    N = L.*L
    α = 2.0

    function get_s_and_h_plots(nr, t)
        files = [load_file(l, α, t, nr) for l in L]
        seeds = zip(L,[f["nr_seeds_used"] for f in files])
        [println("$nr, $t: ($L, $s)") for (L,s) in seeds]
        s = [f["average_spanning_cluster_size"] for f in files]
        h = [f["average_spanning_cluster_perimiter"] for f in files]

        std_s = [f["std_spanning_cluster_size"] for f in files]
        std_h = [f["std_spanning_cluster_perimiter"] for f in files]

        log_L = log2.(L)
        s = log2.(s.± std_s)
        h = log2.(h.± std_s)

        fit_s = get_fit(log_L, Measurements.value.(s))
        fit_h = get_fit(log_L, Measurements.value.(h))
        slope_s = round(fit_s[2], digits=2)
        slope_h = round(fit_h[2], digits=2)

        s = s .± log2.(std_s)
        h = h .± log2.(std_h)

        m_shape = :cross
        default(color=:black, markersize=3)
        s_plot = scatter(log_L, s, markershape=m_shape, label=L"Cluster size ($s$)",
                        legend=:topleft, xlabel=L"log$_2(L)$", ylabel=L"log$_2(s)$",
                        title=latexstring("$nr: \$D_s=$slope_s\$"))
        plot!(log_L, f_lin(log_L, fit_s), label=L"Fit ($D_s$)", color=:black)
        h_plot = scatter(log_L, h, markershape=m_shape, label=L"Perimiter length ($h$)",
                        legend=:topleft, xlabel=L"log$_2(L)$", ylabel=L"log$_2(h)$",
                        title=latexstring("$nr: \$D_h=$slope_h\$"))
        plot!(log_L, f_lin(log_L, fit_h), label=L"Fit ($D_h$)")
        #yticks!(h_plot, yticks(s_plot)[1])
        return s_plot, h_plot
    end

    for t_ in t
        @assert NR[1] == "SNR"
        SNR_s_plot, SNR_h_plot = get_s_and_h_plots(NR[1], t_)
        @assert NR[2] == "UNR"
        UNR_s_plot, UNR_h_plot = get_s_and_h_plots(NR[2], t_)

        l = @layout [
            A B;
            C D
        ]
        plot(SNR_s_plot, SNR_h_plot, UNR_s_plot, UNR_h_plot, size=(800, 600), layout = l,
            plot_title=latexstring("Dimensionality: \$t=$t_\$"))

        savefig("plots/Graphs/dimension_t=$t_.pdf")
    end
end

function get_s_and_h_p_dimension()
    
    L = [8,16,32,64,128,256,512]
    t = (0:9)./10
    NR = ["SNR", "UNR"]
    N = L.*L
    α = 2.0

    function get_s_and_h_p_dimension(nr, t)
        files = [load_file(l, α, t, nr) for l in L]
        #seeds = last(files)["nr_seeds_used"]
        #println(seeds)
        s = [f["average_spanning_cluster_size"] for f in files]
        h = [f["average_spanning_cluster_perimiter"] for f in files]

        std_s = [f["std_spanning_cluster_size"] for f in files]
        std_h = [f["std_spanning_cluster_perimiter"] for f in files]

        log_L = log2.(L)
        s = log2.(s.± std_s)
        h = log2.(h.± std_s)

        fit_s = get_fit(log_L, Measurements.value.(s))
        fit_h = get_fit(log_L, Measurements.value.(h))
        slope_s = round(fit_s[2], digits=2)
        slope_h = round(fit_h[2], digits=2)

        s = s .± log2.(std_s)
        h = h .± log2.(std_h)

        return slope_s, slope_h
    end

    SNR_s_slope = zeros(Float64, length(t))
    SNR_h_slope = zeros(Float64, length(t))
    UNR_s_slope = zeros(Float64, length(t))
    UNR_h_slope = zeros(Float64, length(t))
    for i in eachindex(t)
        SNR_s_slope[i], SNR_h_slope[i] = get_s_and_h_p_dimension(NR[1], t[i])
        UNR_s_slope[i], UNR_h_slope[i] = get_s_and_h_p_dimension(NR[2], t[i])
    end

    plot([t,t,t,t], [SNR_s_slope, SNR_h_slope, UNR_s_slope, UNR_h_slope], size=(400, 300),
        plot_title="Dimensionality", labels=[L"SNR $D_s$" L"SNR $D_h$" L"UNR $D_s$" L"UNR $D_h$"],
        legend=:right, xlabel=L"t", ylabel="Dimensionality", linestyle=[:solid :solid :dash :dash])

    savefig("plots/Graphs/dimension.pdf")
end

function get_gyration_radii(L, nr, t, files)
    bundles = [get_fb(L, nr=nr, t=t) for _ in eachindex(files)]
    


function plot_dimensions_over_t_with_radius_of_gyration()
    
    L = [8,16,32,64,128,256,512]
    t = (0:9)./10
    NR = ["SNR", "UNR"]
    N = L.*L
    α = 2.0

    function get_s_and_gyration(nr, t)
        files = [load_file(l, α, t, nr) for l in L]

        s = [f["average_spanning_cluster_size"] for f in files]
        bundles = [f["average_spanning_cluster_perimiter"] for f in files]

        std_s = [f["std_spanning_cluster_size"] for f in files]
        std_h = [f["std_spanning_cluster_perimiter"] for f in files]

        log_L = log2.(L)
        s = log2.(s.± std_s)
        h = log2.(h.± std_s)

        fit_s = get_fit(log_L, Measurements.value.(s))
        fit_h = get_fit(log_L, Measurements.value.(h))
        slope_s = round(fit_s[2], digits=2)
        slope_h = round(fit_h[2], digits=2)

        s = s .± log2.(std_s)
        h = h .± log2.(std_h)

        return slope_s, slope_h
    end

    SNR_s_slope = zeros(Float64, length(t))
    SNR_h_slope = zeros(Float64, length(t))
    UNR_s_slope = zeros(Float64, length(t))
    UNR_h_slope = zeros(Float64, length(t))
    for i in eachindex(t)
        SNR_s_slope[i], SNR_h_slope[i] = get_s_and_h_p_dimension(NR[1], t[i])
        UNR_s_slope[i], UNR_h_slope[i] = get_s_and_h_p_dimension(NR[2], t[i])
    end

    plot([t,t,t,t], [SNR_s_slope, SNR_h_slope, UNR_s_slope, UNR_h_slope], size=(400, 300),
        plot_title="Dimensionality", labels=[L"SNR $D_s$" L"SNR $D_h$" L"UNR $D_s$" L"UNR $D_h$"],
        legend=:right, xlabel=L"t", ylabel="Dimensionality", linestyle=[:solid :solid :dash :dash])

    savefig("plots/Graphs/dimension_using_radius_of_gyration.pdf")
end

plot_dimensions_over_t()
plot_dimensions_over_t_with_radius_of_gyration()
plot_dimension_thing()



println("Saved plot!")