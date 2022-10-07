using Plots
using JLD2
using LaTeXStrings

include("../support/ploting_settings.jl")

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*".jld2"

global_path = "data/"

distributions(t) = [ "t=$t Uniform", "t=$t Uniform SNR", "t=$t Uniform CNR"]
desired_data = [
    "average_nr_clusters",
    "average_largest_cluster",
    "average_largest_perimiter",
    "average_most_stressed_fiber",
    "nr_seeds_used",
]

function file(global_path, L, distribution)
    # We only extract the data we want
    fileDict = Dict()
    jldopen(full_name(global_path, L, distribution), "r") do f
        for data_key in desired_data
            fileDict[data_key] = f[data_key]
        end
    end
    return fileDict
end

L = 128
N = L.*L
k_N = [1:n for n in N]./N

t = 1:9
lables = permutedims(["UNR", "SNR", "CNR"])
files = [[file(global_path, L, distribution) for distribution in distributions(t_)] for t_ in t ./10]
seeds = files[1][1]["nr_seeds_used"]
#Files [t] [distribution] [value]

max_number_of_clusters = [[1/(maximum(files[t_][i]["average_nr_clusters"])*N) for t_ in t] for i in 1:length(distributions(1))]
k_where_nr_clusters_is_max = [[argmax(files[t_][i]["average_nr_clusters"]) / N for t_ in t] for i in 1:length(distributions(1))]
sigma_c = [[maximum(files[t_][i]["average_most_stressed_fiber"]) for t_ in t][1:4] for i in 1:length(distributions(1))]
k_where_sigma_c_is_max = [[argmax(files[t_][i]["average_most_stressed_fiber"]) / N for t_ in t][1:4] for i in 1:length(distributions(1))]
largest_cluster_size = [[files[t_][i]["average_largest_cluster"][floor(Int64, N/3)] for t_ in t] for i in 1:length(distributions(1))]
largest_perimeter = [[files[t_][i]["average_largest_perimiter"][floor(Int64, N/3)] for t_ in t] for i in 1:length(distributions(1))]


default(markersize=3)

max_number_of_clusters_plot = plot(t ./ 10, max_number_of_clusters , label = lables, legend=:topleft, marker=:circle,
                    xlabel=L"t_0", ylabel=L"1/N^{\mathrm{max}}_c", title="Maximum of number of clusters")

k_where_nr_clusters_is_max_plot = plot(t  ./ 10, k_where_nr_clusters_is_max, label = lables, marker=:circle,
                    xlabel=L"t_0", ylabel=L"k/N"*" at "*L"N^{\mathrm{max}}_c)", title=L"k/N"*" where maximum was reached")
                    
sigma_c_plot = plot(t[1:4]  ./ 10, sigma_c, label = lables, legend=:topleft, marker=:circle,
                    xlabel=L"t_0", ylabel=L"σ_c", title="Maximum of "*L"σ_c")

#k_where_sigma_c_is_max_plot = plot(t[1:4]  ./ 10, k_where_sigma_c_is_max, label = lables,
#                    xlabel=L"t_0", ylabel=L"k/N"*" at "*L"σ_c", title=L"k/N"*" where maximum was reached")


cluster_over_perimiter_size_plot = plot(largest_perimeter, largest_cluster_size, label = lables, legend=:bottomright, marker=:circle, yaxis=:log, xaxis=:log, 
                    xlabel=L"H_{\mathrm{max}}/N", ylabel=L"S_{\mathrm{max}}/N", title="Cluster over perimiter size at "*L"k/N=\frac{1}{3}")
n = 2
cur_colors = theme_palette(:auto)
scatter!(largest_perimeter[1][1:n], largest_cluster_size[1][1:n], color=cur_colors[1], series_annotations = text.(L"t_0=".*latexstring.(t[1:n] ./10).*" ", :right, :bottom, 7), primary=false)
scatter!(largest_perimeter[2][1:n], largest_cluster_size[2][1:n], color=cur_colors[2], series_annotations = text.(" "*L"t_0=".*latexstring.(t[1:n] ./10).*" ", :left, :bottom, 7), primary=false)
scatter!(largest_perimeter[3][1:n], largest_cluster_size[3][1:n], color=cur_colors[3], series_annotations = text.(" "*L"t_0=".*latexstring.(t[1:n] ./10).*" ", :right, :top, 7), primary=false)


lagrest_clusetr_size_plot = plot(t  ./ 10, largest_cluster_size, label = lables, legend=:bottomright, marker=:circle,
                    xlabel=L"t_0", ylabel=L"S_{\mathrm{max}}/N", title="Largest cluster size at "*L"k/N=\frac{1}{3}")

largest_perimeiter_plot = plot(t  ./ 10, largest_perimeter, label = lables, legend=:right, marker=:circle,
                    xlabel=L"t_0", ylabel=L"H_{\mathrm{max}}/N", title="Largest perimiter size at "*L"k/N=\frac{1}{3}")


l = @layout [
    A B; C D; E F
]
plot(max_number_of_clusters_plot, k_where_nr_clusters_is_max_plot, sigma_c_plot, cluster_over_perimiter_size_plot, lagrest_clusetr_size_plot, largest_perimeiter_plot,
    size=(800, 800), layout = l, plot_title=latexstring("Uniform distribution with NHR over regiemes, \$L=$L\$, $seeds samples"), plot_titlevspan=0.1)

savefig("plots/Graphs/Uniform with Neighbourhood rules over different regiemes.pdf")

println("Saved plot!")