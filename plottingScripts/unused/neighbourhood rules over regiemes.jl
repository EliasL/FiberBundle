using Plots
using JLD2
using LaTeXStrings
using Measures

include("ploting_settings.jl")
include("../support/dataManager.jl")

global_path = "data/"

L = 32
N = L*L
α = 2.0
t = vcat((0:20) ./ 50, (5:9) ./ 10)
dists = ["LLS", "CLS"]
lables = permutedims(dists)
files = [[load_file(L, α, t_, nr) for nr in dists] for t_ in t]
seeds = files[1][1]["nr_seeds_used"]
#Files [t] [distribution] [value]

function get_data_from_files(key, f::Function = x -> x; t_index=:, d_index=:)
    try
        return [[f(files[j][i][key][d_index]) for j in 1:length(t)][t_index] for i in 1:length(dists)]
    catch
        return [[f(files[j][i][key]) for j in 1:length(t)][t_index] for i in 1:length(dists)]
    end
end

max_number_of_clusters = get_data_from_files("average_nr_clusters", x -> 1/(maximum(x)*N))
k_where_nr_clusters_is_max = get_data_from_files("average_nr_clusters", x -> argmax(x)/N)
sigma_c_index = 1:25
sigma_c = get_data_from_files("average_most_stressed_fiber", maximum, t_index=sigma_c_index)
k_where_sigma_c_is_max = get_data_from_files("average_most_stressed_fiber", x -> argmax(x)/N)
largest_cluster_size = get_data_from_files("average_largest_cluster", d_index = floor(Int64, N/5))
largest_perimeter = get_data_from_files("average_largest_perimiter", d_index = floor(Int64, N/5))
spanning_cluster_size = get_data_from_files("average_spanning_cluster_size", x -> x / N)
spanning_perimeter = get_data_from_files("average_spanning_cluster_perimiter", x -> x / N)

max_number_of_clusters_plot = plot(t, max_number_of_clusters , label = lables, legend=:bottomright, marker=:circle,
                    xlabel=L"t_0", ylabel=L"1/N^{\mathrm{max}}_c", title="A: Maximum of number of clusters")

k_where_nr_clusters_is_max_plot = plot(t, k_where_nr_clusters_is_max, label = lables, marker=:circle,
                    xlabel=L"t_0", ylabel=L"k/N"*" at "*L"N^{\mathrm{max}}_c", title=L"B: k/N"*" where maximum was reached")
              
sigma_c_plot = plot(t[sigma_c_index], sigma_c, label = lables, legend=:bottomright, marker=:circle, ylims=(0, Inf),
                    xlabel=L"t_0", ylabel=L"σ_c", title="C: Maximum of "*L"σ_c")

#k_where_sigma_c_is_max_plot = plot(t[1:4]  ./ 10, k_where_sigma_c_is_max, label = lables,
#                    xlabel=L"t_0", ylabel=L"k/N"*" at "*L"σ_c", title=L"k/N"*" where maximum was reached")


cluster_over_perimiter_size_plot = plot(spanning_perimeter, spanning_cluster_size, label = lables, legend=:topright,
                    marker=:circle, yaxis=:log, xaxis=:log, 
                    xlabel=L"H_{\mathrm{max}}/N", ylabel=L"S_{\mathrm{max}}/N", title="A: Spanning cluster over perimiter")

r1 = [1,5,9]
r2 = [1,3,5,11]
cur_colors = theme_palette(:auto)
#scatter!(spanning_perimeter[1][s:s+n], spanning_cluster_size[1][s:s+n], color=cur_colors[1], series_annotations = text.(L"t_0=".*latexstring.(t[s:s+n] ./10).*" ", :right, :bottom, 7), primary=false)
scatter!(spanning_perimeter[2][r1], spanning_cluster_size[2][r1], color=cur_colors[2],
    series_annotations = Plots.text.(""*L"t_0=".*latexstring.(t[r1]), :left, :bottom, 7), primary=false)


lagrest_clusetr_size_plot = plot(t, largest_cluster_size, label = lables, legend=:bottomright, marker=:circle, ylims=(0, Inf),
                    xlabel=L"t_0", ylabel=L"S_{\mathrm{max}}/N", title="D: Largest cluster size at "*L"k/N=1/5")

largest_perimeter_plot = plot(t, largest_perimeter, label = lables, legend=:right, marker=:circle, ylims=(0, Inf),
                    xlabel=L"t_0", ylabel=L"H_{\mathrm{max}}/N", title="D: Largest perimiter size at "*L"k/N=1/5")

spanning_cluster_size_plot = plot(t, spanning_cluster_size, label = lables, legend=:topleft, marker=:circle, ylims=(0, Inf),
                    xlabel=L"t_0", ylabel=L"S_{\mathrm{span}}/N", title="B: Spanning cluster size")

spanning_perimeter_plot = plot(t, spanning_perimeter, label = lables, legend=:right, marker=:circle, ylims=(0, Inf),
                    xlabel=L"t_0", ylabel=L"H_{\mathrm{span}}/N", title="C: Spanning perimiter size")

function plot_all()
    l = @layout [
        A B; C D; E F ; G H
    ]
    plot(max_number_of_clusters_plot, k_where_nr_clusters_is_max_plot, sigma_c_plot, cluster_over_perimiter_size_plot, spanning_cluster_size_plot, spanning_perimeter_plot, lagrest_clusetr_size_plot, largest_perimeter_plot,
        size=(1000, 1000), layout = l,
        plot_title=latexstring("Uniform distribution over regiemes, \$L=$L\$, $seeds samples"), plot_titlevspan=0.1)

    savefig("plots/Graphs/Uniform with Neighbourhood rules over different regiemes.pdf")

    println("Saved plot!")
end

function plot_similar()
    l = @layout [
        A; B; C; D;
    ]
    plot(max_number_of_clusters_plot, k_where_nr_clusters_is_max_plot, sigma_c_plot,
        size=(500, 900), layout = l, left_margin=10Plots.mm)#, bottom_margin=-3Plots.mm, right_margin=3Plots.mm)

    savefig("plots/Graphs/NR_similarities.pdf")
end

function plot_differences()
    l = @layout [
        A; B; C; D;
    ]
    plot(cluster_over_perimiter_size_plot, spanning_cluster_size_plot, spanning_perimeter_plot,
    size=(500, 900), layout = l, left_margin=10Plots.mm)#, bottom_margin=-3Plots.mm, right_margin=3Plots.mm)

    savefig("plots/Graphs/NR_differences.pdf")
end


plot_similar()
plot_differences()
plot_all()
println("Saved plots!")