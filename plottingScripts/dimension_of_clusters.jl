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
    NR = ["CLS", "LLS"]
    N = L.*L
    α = 2.0

    function get_s_and_r_plots(nr, t)
        files = [load_file(l, α, t, nr) for l in L]
        seeds = zip(L,[f["nr_seeds_used"] for f in files])
        f = load_file(L[1], α, t, nr, average=false)
        seeds = 0:4000-1
        println("$nr, $t")
        for j in [40, 400, 2000, 4000-1]
            println(std([f["spanning_cluster_perimiter/$i"] for i in 0:j]))
        end
        spanning_cluster_size = [f["spanning_cluster_size/$i"] for i in seeds]
        display(plot(seeds, spanning_cluster_size))
        s = [f["average_spanning_cluster_size"] for f in files]
        h = [f["average_spanning_cluster_perimiter"] for f in files]

        std_s = [f["std_spanning_cluster_size"] for f in files]
        std_h = [f["std_spanning_cluster_perimiter"] for f in files]

        log_L = log2.(L)
        s = log2.(s.± std_s)
        h = log2.(h.± std_h)

        fit_s = get_fit(log_L, Measurements.value.(s))
        fit_h = get_fit(log_L, Measurements.value.(h))
        slope_s = round(fit_s[2], digits=2)
        slope_h = round(fit_h[2], digits=2)

        #s = s .± log2.(std_s)
        #h = h .± log2.(std_h)

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
        @assert NR[1] == "CLS"
        CLS_s_plot, CLS_r_plot = get_s_and_r_plots(NR[1], t_)
        @assert NR[2] == "LLS"
        LLS_s_plot, LLS_r_plot = get_s_and_r_plots(NR[2], t_)

        l = @layout [
            A B;
            C D
        ]
        plot(CLS_s_plot, CLS_r_plot, LLS_s_plot, LLS_r_plot, size=(800, 600), layout = l,
            plot_title=latexstring("Dimensionality: \$t=$t_\$"))

        savefig("plots/Graphs/dimension_t=$t_.pdf")
    end
end

function plot_dimensions_over_t()
    
    L = [8,16,32,64,128,256,512]
    t = (0:9)./10
    NR = ["CLS", "LLS"]
    N = L.*L
    α = 2.0

    function get_s_and_r_p_dimension(nr, t)
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

        #s = s .± log2.(std_s)
        #h = h .± log2.(std_h)

        return slope_s, slope_h
    end

    CLS_s_slope = zeros(Float64, length(t))
    CLS_h_slope = zeros(Float64, length(t))
    LLS_s_slope = zeros(Float64, length(t))
    LLS_h_slope = zeros(Float64, length(t))
    for i in eachindex(t)
        CLS_s_slope[i], CLS_h_slope[i] = get_s_and_r_p_dimension(NR[1], t[i])
        LLS_s_slope[i], LLS_h_slope[i] = get_s_and_r_p_dimension(NR[2], t[i])
    end

    plot([t,t,t,t], [CLS_s_slope, CLS_h_slope, LLS_s_slope, LLS_h_slope], size=(400, 300),
        plot_title="Dimensionality", labels=[L"CLS $D_s$" L"CLS $D_h$" L"LLS $D_s$" L"LLS $D_h$"],
        legend=:right, xlabel=L"t", ylabel="Dimensionality", linestyle=[:solid :solid :dash :dash])

    savefig("plots/Graphs/dimension.pdf")
end

function get_gyration_radii(L, nr, t, bulk_files)

    try
        f = load("data/gyration_data/")
        CLS_r_slope = f["CLS_r_slope"]
        LLS_r_slope = f["LLS_r_slope"]
        CLS_s_slope = f["CLS_s_slope"]
        LLS_s_slope = f["LLS_s_slope"]
    catch
        for i in eachindex(t)
            CLS_s_slope[i], CLS_r_slope[i] = get_s_and_gyration(NR[1], t[i])
            LLS_s_slope[i], LLS_r_slope[i] = get_s_and_gyration(NR[2], t[i])
        end
        jldopen("data/gyration_data.jld2", "w") do file
            file["CLS_r_slope"] = CLS_r_slope
            file["LLS_r_slope"] = LLS_r_slope
            file["CLS_s_slope"] = CLS_s_slope
            file["LLS_s_slope"] = LLS_s_slope
        end
    end

    bundles = [get_fb(L[i], nr=nr, t=t, without_storage=true) for i in eachindex(bulk_files)]
    r = []
    std_r = []
    for (b::FB, file) in zip(bundles, bulk_files)
        temp_r = []
        for seed in file["seeds_used"]
            #println(view(file["break_sequence/$seed"],1:file["last_step/$seed"]))
            
            break_fiber_list!(view(file["break_sequence/$seed"],1:file["last_step/$seed"]), b)
            update_σ!(b)
            #TODO ASSUMPTION MAXIMUM CLUSTER IS SPANNING CLUSTER
            single_r = maximum(find_radius_of_gyration(b))
            push!(temp_r, single_r)
            healBundle!(b)
        end
        push!(r, mean(temp_r))
        push!(std_r, std(temp_r))
    end
    return r, std_r
end


function plot_dimensions_over_t_with_radius_of_gyration()
    
    L = [8,16,32,64,128,256,512]
    t = (0:9)./10
    NR = ["CLS", "LLS"]
    N = L.*L
    α = 2.0

    function get_s_and_gyration(nr, t)
        files = [load_file(l, α, t, nr) for l in L]
        bulk_files = [load_file(l, α, t, nr,average=false) for l in L]

        s = [f["average_spanning_cluster_size"] for f in files]
        std_s = [f["std_spanning_cluster_size"] for f in files]
        r, std_r = get_gyration_radii(L, nr, t, bulk_files)

        log_L = log2.(L)
        s = log2.(s.± std_s)
        r = log2.(r.± std_r)

        fit_s = get_fit(log_L, Measurements.value.(s))
        fit_r = get_fit(log_L, Measurements.value.(r))
        slope_s = round(fit_s[2], digits=2)
        slope_r = round(fit_r[2], digits=2)

        return slope_s, slope_r
    end

    CLS_s_slope = zeros(Float64, length(t))
    CLS_r_slope = zeros(Float64, length(t))
    LLS_s_slope = zeros(Float64, length(t))
    LLS_r_slope = zeros(Float64, length(t))

    for i in eachindex(t)
        CLS_s_slope[i], CLS_r_slope[i] = get_s_and_gyration(NR[1], t[i])
        LLS_s_slope[i], LLS_r_slope[i] = get_s_and_gyration(NR[2], t[i])
    end

    plot([t,t,t,t], [CLS_s_slope, CLS_r_slope, LLS_s_slope, LLS_r_slope], size=(400, 300),
        plot_title="Dimensionality", labels=[L"CLS $D_s$" L"CLS $D_r$" L"LLS $D_s$" L"LLS $D_r$"],
        legend=:right, xlabel=L"t", ylabel="Dimensionality", linestyle=[:solid :solid :dash :dash])

    savefig("plots/Graphs/dimension_using_radius_of_gyration.pdf")
end

function plot_dimensions_with_radius_of_gyration()
    
    L = [8,16,32,64,128,256,512]
    t = (0:9)./10
    NR = ["CLS", "LLS"]
    N = L.*L
    α = 2.0

    function get_s_and_gyration_plot(nr, t)
        files = [load_file(l, α, t, nr) for l in L]
        bulk_files = [load_file(l, α, t, nr,average=false) for l in L]

        s = [f["average_spanning_cluster_size"] for f in files]
        std_s = [f["std_spanning_cluster_size"] for f in files]
        r, std_r = get_gyration_radii(L, nr, t, bulk_files)

        log_L = log2.(L)
        s = log2.(s.± std_s)
        r = log2.(r.± std_r)

        fit_s = get_fit(log_L, Measurements.value.(s))
        fit_r = get_fit(log_L, Measurements.value.(r))
        slope_s = round(fit_s[2], digits=2)
        slope_r = round(fit_r[2], digits=2)

        
        m_shape = :cross
        default(color=:black, markersize=3)
        s_plot = scatter(log_L, s, markershape=m_shape, label=L"Cluster size ($s$)",
                        legend=:topleft, xlabel=L"log$_2(L)$", ylabel=L"log$_2(s)$",
                        title=latexstring("$nr: \$D_s=$slope_s\$"))
        plot!(log_L, f_lin(log_L, fit_s), label=L"Fit ($D_s$)", color=:black)
        r_plot = scatter(log_L, r, markershape=m_shape, label=L"Perimiter length ($r$)",
                        legend=:topleft, xlabel=L"log$_2(L)$", ylabel=L"log$_2(r)$",
                        title=latexstring("$nr: \$D_r=$slope_r\$"))
        plot!(log_L, f_lin(log_L, fit_r), label=L"Fit ($D_r$)")
        #yticks!(h_plot, yticks(s_plot)[1])
        return s_plot, h_plot
    end

    CLS_s_slope = zeros(Float64, length(t))
    CLS_r_slope = zeros(Float64, length(t))
    LLS_s_slope = zeros(Float64, length(t))
    LLS_r_slope = zeros(Float64, length(t))


    for t_ in t
        @assert NR[1] == "CLS"
        CLS_s_plot, CLS_r_plot = get_s_and_r_plots(NR[1], t_)
        @assert NR[2] == "LLS"
        LLS_s_plot, LLS_r_plot = get_s_and_r_plots(NR[2], t_)

        l = @layout [
            A B;
            C D
        ]
        plot(CLS_s_plot, CLS_r_plot, LLS_s_plot, LLS_r_plot, size=(800, 600), layout = l,
            plot_title=latexstring("Dimensionality: \$t=$t_\$"))

        savefig("plots/Graphs/dimension_t=$t_.pdf")
    end


    plot([t,t,t,t], [CLS_s_slope, CLS_r_slope, LLS_s_slope, LLS_r_slope], size=(400, 300),
        plot_title="Dimensionality", labels=[L"CLS $D_s$" L"CLS $D_r$" L"LLS $D_s$" L"LLS $D_r$"],
        legend=:right, xlabel=L"t", ylabel="Dimensionality", linestyle=[:solid :solid :dash :dash])

    savefig("plots/Graphs/dimension_using_radius_of_gyration.pdf")
end

plot_dimensions_over_t()
plot_dimensions_over_t_with_radius_of_gyration()
plot_dimension_thing()



println("Saved plot!")