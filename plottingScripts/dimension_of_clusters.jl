using Plots, Measurements
using JLD2
using LaTeXStrings
using LsqFit
using ProgressMeter


include("ploting_settings.jl")
include("../support/dataManager.jl")
include("../support/bundleAnalasys.jl")
include("showBundle.jl")

f_lin(x,p) = p[1] .+ Measurements.value.(x) .* p[2]

function get_fit(x, y)
    p0_lin = [0.0, 0.0]
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
    i=('s' in label) ? 1 : 2
    c = theme_palette(:auto)[i]

    scatter!(x, y, markershape=marker, label=label, linewidth=1.2,
    markerstrokecolor=c,legend=:topleft, xlabel=x_label,
    markersize=5, markerstrokewidth=1.2)
    #   y fit
    plot!(Measurements.value.(x), f_lin(x, fit), label="", color=c,
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
            push!(labels, latexstring("D_s^{$nr}"))
            h_plot, h_slope = get_plot_and_slope!(p, r, h, " "*L"h",
            :dot, :dtriangle, x_label=L"log$_2(R)$")
            push!(plots, h_plot)
            slopes[j, i*2] =  h_slope
            push!(labels, latexstring("D_h^{$nr}"))
            s_rounded = round(s_slope, digits=2)
            h_rounded = round(h_slope, digits=2)
            if i==1
                ylabel!(L"log$_2(s)$, log$_2(h)$")
            end
            title!(latexstring("D_s^{$nr}=")*"$s_rounded "*
                   latexstring("D_h^{$nr}=")*"$h_rounded")
            end
        
         # Plot for each t
        plot(plots..., size=(400, 200), layout = (length(NR)),)
        #plot_title=latexstring("Dimensionality: \$t_0=$(t)\$"))    
        savefig("plots/Graphs/SingleT/new_dist_dimension_t=$(t)_L=$(L[1])-$(L[end]).pdf")

    end
    #scatter!(x_newData/N, y, markerstrokecolor=colors, markercolor=:transparent, label=nothing, markershape=:diamond, markersize=5, markerstrokewidth=1)
    scatter(ts, slopes, labels=permutedims(labels), markerstrokecolor=permutedims(theme_palette(:auto)[1:4]), markercolor=:transparent, markersize=5, markershape=[:utriangle :dtriangle :diamond :star4 :star6 :pentagon],
        size=(300, 250), legend=:left, xlabel=L"t_0", ylabel=L"D")
    #println(slopes)
    savefig("plots/Graphs/new_dist_dimension_L=$(L[1])-$(L[end]).pdf")
    
end

function plot_dimension_box_thing(l, ts, α)
    NR = ["CLS", "LLS"]
    ms = [:utriangle, :dtriangle, :star4, :diamond]
    slopes = zeros(Float64, (length(ts) , length(NR), 2))
    c = permutedims(theme_palette(:auto)[1:4])
    for (j,t) in enumerate(ts)
        
        # Get plot and slope
        p = plot(size=(320, 320), ylabel=L"\log_2(s), \log_2(h)",
                xlabel=L"\log_2(\epsilon)")
        for (i, nr) in enumerate(NR)
            # we don't like the last data points
            end_datapoint = 6 #Max 9, 9=256

            start_datapoint = 2 #Min 1
            count = get_box_count(l, nr, t, α)
            s_count, h_count = log2.(count[1][start_datapoint:end_datapoint]), log2.(count[2][start_datapoint:end_datapoint])
            nr_sizes = round(Int64, log2(l))-1
            box_sizes = start_datapoint-1:end_datapoint-1 #log2(box_sizes) actually
            
            s_fit = get_fit(box_sizes, s_count)
            s_slope = s_fit[2]
            s_rounded_slope = round(s_slope, digits=2)
            slopes[j, i, 1] = -s_slope
            scatter!(box_sizes, s_count, label=latexstring("D_s^{$nr}")*": slope $s_rounded_slope", markershape=ms[i*2-1],
            markerstrokecolor=c[i*2-1])
            plot!(box_sizes, f_lin(box_sizes, s_fit), label="", c=c[i*2-1])

            h_fit = get_fit(box_sizes, h_count)
            h_slope = h_fit[2]
            h_rounded_slope = round(h_slope, digits=2)
            slopes[j, i, 2] = -h_slope
            scatter!(box_sizes, h_count, label="$nr "*latexstring("D_h^{$nr}")*": slope $h_rounded_slope", markershape=ms[i*2],
                markerstrokecolor=c[i*2])
            plot!(box_sizes, f_lin(box_sizes, h_fit), label="", c=c[i*2])
         # Plot for each t
        end
        #plot_title=latexstring("Dimensionality: \$t_0=$(t)\$"))    
        savefig("plots/Graphs/SingleT/box_dimension_t=$(t)_L=$l.pdf")
    end
    #scatter!(x_newData/N, y, markerstrokecolor=colors, markercolor=:transparent, label=nothing, markershape=:diamond, markersize=5, markerstrokewidth=1)
    s = hcat(slopes[:, 1, :], slopes[:, 2, :])
    scatter(ts, s,
        size=(300, 250), legend=:left, xlabel=L"t_0", ylabel=L"D", markershape=[:utriangle :dtriangle :diamond :star4],
        markerstrokecolor=c, label=[L"$D_s^{CLS}$" L"$D_h^{CLS}$" L"$D_s^{LLS}$" L"$D_h^{LLS}$"])
    #println(slopes)
    savefig("plots/Graphs/box_dimension_L=$l.pdf")
    
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

function do_box_count(l, nr, t, file)
    seeds = file["seeds_used"]
    nr_seeds = length(seeds)
    s_counts = zeros(round(Int64, log2(l)))
    h_counts = zeros(round(Int64, log2(l)))
    for seed in seeds
        b = get_bundle_from_file(file, l, nr=nr, t=t, spanning=true, seed=seed)

        counts = box_counting(b)
        s_counts .+= counts[1]
        h_counts .+= counts[2]
    end
    return s_counts ./ nr_seeds, h_counts ./ nr_seeds
end

function get_box_count(l, nr, t, α)
    if !isdir("newData/box_counting_data/")
        mkpath("newData/box_counting_data/")
    end
    file_name = get_file_name(l, α, t, nr, dist, data_path=data_path)
    if isfile(file_name)
        if isfile("newData/box_counting_data/$(l)_$(nr)_$(t)_boxes.jld2")
            f = load("newData/box_counting_data/$(l)_$(nr)_$(t)_boxes.jld2")
            count = f["$(l)_$(nr)_$(t)_boxes"]
        else
            bulk_file = load_file(l, α, t, nr, dist, average=false, data_path=data_path)
            count = do_box_count(l, nr, t, bulk_file)
            # Save data
            jldopen("newData/box_counting_data/$(l)_$(nr)_$(t)_boxes.jld2", "w") do file
                file["$(l)_$(nr)_$(t)_boxes"] = count
            end
        end
    else 
        @error "No file found! $file_name" 
    end
    return count
end

function get_gyration_radii(L, nr, t, α)
    R = []
    S = []
    H = []
    std_R = []
    std_S = []
    std_H = []
    if !isdir("newData/gyration_data/")
        mkpath("newData/gyration_data/")
    end
    for l in L
        file_name = get_file_name(l, α, t, nr, dist, data_path=data_path)
        if isfile(file_name)
            if isfile("newData/gyration_data/$(l)_$(nr)_$(t)_r.jld2")
                f = load("newData/gyration_data/$(l)_$(nr)_$(t)_r.jld2")
                r, s, h = f["$(l)_$(nr)_$(t)_r"]
            else
                bulk_file = load_file(l, α, t, nr,dist, average=false, data_path=data_path)
                r, s, h = calculate_gyration_radi(l, nr, t, bulk_file)
                # Save data
                jldopen("newData/gyration_data/$(l)_$(nr)_$(t)_r.jld2", "w") do file
                    file["$(l)_$(nr)_$(t)_r"] = (r, s, h)
                end
            end
            push!(R, mean(r))
            push!(std_R, std(r))
            push!(S, mean(s))
            push!(std_S, std(s))
            push!(H, mean(h))
            push!(std_H, std(h))
            #println("$l, $t, $(length(r))")
        else 
            @error "No file found! $file_name" 
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

L = [16, 32, 64, 128, 256, 512]#512, 1024]
α = 2.0
dist = "ConstantAverageUniform"
data_path = "newData/"
#t = vcat((1:9) ./ 10)
#t = vcat((0:1) ./ 10, (10:20) ./ 50, (5:9) ./ 10)
#t = vcat((0:20) ./ 50, (5:9) ./ 10)
#t = vcat(0.05:0.05:0.20, 0.25:0.01:0.5, [0.6, 0.7, 0.8, 0.1])
t = vcat(0.05:0.05:0.20, 0.25:0.01:0.5)
plot_dimension_thing(L, t, α)
plot_dimension_box_thing(512, t, α)


println("Saved plot!")