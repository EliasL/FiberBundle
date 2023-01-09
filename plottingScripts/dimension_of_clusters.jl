using Plots, Measurements
using JLD2
using LaTeXStrings
using LsqFit
using ProgressMeter


include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/inertia.jl")
include("showBundle.jl")

f_lin(x,p) = p[1] .+ Measurements.value.(x) .* p[2]

function get_fit(x, y)
    p0_lin = [0.0, 1.0]
    return curve_fit(f_lin, x, y, p0_lin).param
end


function get_plot_and_slope!(p, x, y, label, linestyle, marker; x_label=L"log$_2(L)$")
    # Converting to log
    x = log2.(x)
    y = log2.(y)

    # Fit to line
    fit = get_fit(Measurements.value.(x), Measurements.value.(y))
    slope = fit[2]
    rounded_slope = round(slope, digits=2)
    
    #= # Get slope error
    slope_err_x, slope_err_y, max_slope, min_slope = uncertainty_in_slope(x, y)
    slope_err = round(abs(slope-max_slope), digits=2)

    # Plot
    #   Slope max_min_error
    p = plot(slope_err_x, slope_err_y, label=nothing, color=:red,
    linewidth=1, linestyle=:dash, markersize=3, markershape=:none)
    #   y values =#
    
    scatter!(x, y, markershape=marker, color=:black, label=label, linewidth=1,
    markerstrokecolor=:black,legend=:topleft, xlabel=x_label,
    markersize=3, markerstrokewidth=1)
    #   y fit
    plot!(Measurements.value.(x), f_lin(x, fit), label="", color=:black,
    linewidth=1, linestyle=linestyle)

    return p, slope# ± slope_err
end

function plot_dimension_thing(L, ts, α)
    NR = ["CLS", "LLS"]
    slopes = zeros(Float64, (length(ts) ,2*length(NR)))
    labels = []


    for (j,t) in enumerate(ts)
        
        # Get plot and slope
        plots = []

        for (i, nr) in enumerate(NR)
            p = plot(titlefontsize=9, xlims=(1,Inf), ylims=(-Inf, 18))
            r, s, h, std_r, std_s, std_h = get_gyration_radii(L, nr, t, α)

            #= s = s .± std_s
            h = h .± std_h
            r = r .± std_r =#
            s_plot, s_slope = get_plot_and_slope!(p, r, s, " "*L"s",
            :solid, :diamond, x_label=L"log$_2(R)$")
            slopes[j, i*2-1] = s_slope
            push!(labels, "$nr "*L"D_s")
            h_plot, h_slope = get_plot_and_slope!(p, r, h, " "*L"h",
            :dot, :dtriangle, x_label=L"log$_2(R)$")
            push!(plots, h_plot)
            slopes[j, i*2] =  h_slope
            push!(labels, "$nr "*L"D_h")
            s_rounded = round(s_slope, digits=2)
            h_rounded = round(h_slope, digits=2)
            if i==1
                ylabel!(L"log$_2(s)$, log$_2(h)$")
            end
            title!("$nr: "*L"D_s="*"$s_rounded "*L"D_h="*"$h_rounded")
            end
        
        # Plot for each t
        plot(plots..., size=(400, 200), layout = (length(NR)),)
        #plot_title=latexstring("Dimensionality: \$t_0=$(t)\$"), )    
        savefig("plots/Graphs/SingleT/dimension_t=$(t)_L=$(L[1])-$(L[end]).pdf")

    end

    plot(ts, slopes, labels=permutedims(labels), markershape=[:utriangle :dtriangle :star4 :diamond :star6 :pentagon],
        legendfontsize=8, size=(400, 300), 
        legend=:right, xlabel=L"t_0", ylabel=L"D")
    savefig("plots/Graphs/dimension_L=$(L[1])-$(L[end]).pdf")
    
end

function calculate_gyration_radi(l, nr, t, file)
    b = get_fb(l, 0, nr=nr, t=t, without_storage=true)
    r = []
    s = []
    h = []
    println("Calculating radius of gyration $l $nr $t")
    for seed in file["seeds_used"]
        break_fiber_list!(view(file["break_sequence/$seed"],1:file["spanning_cluster_step/$seed"]-1), b)
        update_σ!(b)
        id = argmax(b.cluster_size)
        push!(r, find_radius_of_gyration(b)[id])
        push!(s, b.cluster_size[id])
        push!(h, b.cluster_outline_length[id])
        healBundle!(b)
    end
    return r, s, h
end

function get_gyration_radii(L, nr, t, α)
    R = []
    S = []
    H = []
    std_R = []
    std_S = []
    std_H = []
    if !isdir("data/gyration_data/")
        mkpath("data/gyration_data/")
    end
    for l in L
        if isfile(get_file_name(l, α, t, nr))
            if isfile("data/gyration_data/$(l)_$(nr)_$(t)_r.jld2")
                f = load("data/gyration_data/$(l)_$(nr)_$(t)_r.jld2")
                r, s, h = f["$(l)_$(nr)_$(t)_r"]
            else
                bulk_file = load_file(l, α, t, nr,average=false)
                r, s, h = calculate_gyration_radi(l, nr, t, bulk_file)
                # Save data
                jldopen("data/gyration_data/$(l)_$(nr)_$(t)_r.jld2", "w") do file
                    file["$(l)_$(nr)_$(t)_r"] = (r, s, h)
                end
            end
            push!(R, mean(r))
            push!(std_R, std(r))
            push!(S, mean(s))
            push!(std_S, std(s))
            push!(H, mean(h))
            push!(std_H, std(h))
        end

    end
    return R, S, H, std_R, std_S, std_H
end

function uncertainty_in_slope(x, v)

    # TODO: Replace this with student t stuff
    #https://tma4245.math.ntnu.no/enkel-line%C3%A6r-regresjon/inferens-regresjonsparametrene/?fbclid=IwAR3SuPlGuPJ-kTLwaHR0SVxOe1X9FRLvk9iFs8qIa_tX6oM03HfB_r_AXBU
    x = Measurements.value.(x)
    p1_min = (x[1], v[1].val-v[1].err)
    p1_max = (x[1], v[1].val+v[1].err)
    p2_min = (x[end], v[end].val-v[end].err)
    p2_max = (x[end], v[end].val+v[end].err)
    x = [p1_min[1] p1_max[1]; p2_max[1] p2_min[1]]
    y = [p1_min[2] p1_max[2]; p2_max[2] p2_min[2]]
    #[1 1; 2 2], [1 2; 2 3],

    max_slope = (p2_max[2] - p1_min[2])/(p2_max[1] - p1_min[1])
    min_slope = (p2_min[2] - p1_max[2])/(p2_min[1] - p1_max[1])

    return x, y, max_slope, min_slope
end

L = [8, 16, 32, 64, 128, 256, 512, 1024]
α = 2.0

#t = vcat((1:9) ./ 10)
#t = vcat((0:1) ./ 10, (10:20) ./ 50, (5:9) ./ 10)
t = vcat((0:20) ./ 50, (5:9) ./ 10)
#t = [0.0]
#plot_dimensions_over_t(L, t)
#plot_dimensions_over_t_with_radius_of_gyration(L, t)
#plot_dimensions_with_radius_of_gyration(L, t)
plot_dimension_thing(L, t, α)



println("Saved plot!")