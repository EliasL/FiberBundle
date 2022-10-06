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
lables = permutedims(["Uniform", "Uniform SNR", "Uniform CNR"])
files = [[file(global_path, L, distribution) for distribution in distributions(t_)] for t_ in t ./10]
seeds = files[1][1]["nr_seeds_used"]
#Files [t] [distribution] [value]

max_number_of_clusters = [[1/(maximum(files[t_][i]["average_nr_clusters"])*N) for t_ in t] for i in 1:length(distributions(1))]
k_where_nr_clusters_is_max = [[argmax(files[t_][i]["average_nr_clusters"]) / N for t_ in t] for i in 1:length(distributions(1))]
sigma_c = [[maximum(files[t_][i]["average_most_stressed_fiber"]) for t_ in t][1:4] for i in 1:length(distributions(1))]
k_where_sigma_c_is_max = [[argmax(files[t_][i]["average_most_stressed_fiber"]) / N for t_ in t][1:4] for i in 1:length(distributions(1))]
largest_cluster_size = [[files[t_][i]["average_largest_cluster"][floor(Int64, N/3)] for t_ in t] for i in 1:length(distributions(1))]
largest_perimeiter = [[maximum(files[t_][i]["average_largest_perimiter"]) for t_ in t] for i in 1:length(distributions(1))]

max_number_of_clusters_plot = plot(t ./ 10, max_number_of_clusters , label = lables, legend=:topleft,
                    xlabel=L"t_0", ylabel=L"1/N^{\mathrm{max}}_c", title="Maximum of number of clusters")

k_where_nr_clusters_is_max_plot = plot(t  ./ 10, k_where_nr_clusters_is_max, label = lables,
                    xlabel=L"t_0", ylabel=L"k/N"*" at "*L"N^{\mathrm{max}}_c)", title=L"k/N"*" where maximum was reached")
                    
sigma_c_plot = plot(t[1:4]  ./ 10, sigma_c, label = lables, legend=:topleft,
                    xlabel=L"t_0", ylabel=L"σ_c", title="Maximum of "*L"σ_c")

k_where_sigma_c_is_max_plot = plot(t[1:4]  ./ 10, k_where_sigma_c_is_max, label = lables,
                    xlabel=L"t_0", ylabel=L"k/N"*" at "*L"σ_c", title=L"k/N"*" where maximum was reached")

lagrest_clusetr_size_plot = plot(t  ./ 10, largest_cluster_size, label = lables,
                    xlabel=L"t_0", ylabel=L"S_{\mathrm{max}}/N", title="Largest cluster size at "*L"k/N=\frac{1}{3}")

largest_perimeiter_plot = plot(t  ./ 10, largest_perimeiter, label = lables,
                    xlabel=L"t_0", ylabel=L"H_{\mathrm{max}}/N", title="Largest perimiter size")


l = @layout [
    A B; C D; E F
]
plot(max_number_of_clusters_plot, k_where_nr_clusters_is_max_plot, sigma_c_plot, k_where_sigma_c_is_max_plot, lagrest_clusetr_size_plot, largest_perimeiter_plot,
    size=(800, 800), layout = l, plot_title=latexstring("Cluster size over regiemes, \$L=$L\$, $seeds samples"), plot_titlevspan=0.1)

savefig("plots/Graphs/Uniform with Neighbourhood rules over different regiemes.pdf")

println("Saved plot!")