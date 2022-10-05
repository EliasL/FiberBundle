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

k_where_nr_clusters_is_max = [[argmax(files[t_][i]["average_nr_clusters"]) / N for t_ in t] for i in 1:length(distributions(1))]
max_number_of_clusters = [[1/(maximum(files[t_][i]["average_nr_clusters"])*N) for t_ in t] for i in 1:length(distributions(1))]
lagrest clusetr size
largest perimeiter
sigma
time of sigma_c


nr_clusters_plot = plot(t ./ 10, max_number_of_clusters , label = lables, legend=:topleft,
                    xlabel=L"t_0", ylabel=L"1/N^{\mathrm{max}}_c", title="Maximum of number of clusters")

largest_cluster_plot = plot(t  ./ 10, k_where_nr_clusters_is_max, label = lables,
                    xlabel=L"t_0", ylabel=L"k/N"*" at "*L"N^{\mathrm{max}}_c)", title=L"k/N"*" where maximum was reached")

plot(nr_clusters_plot, largest_cluster_plot, 
    plot_title=latexstring("Cluster size over regiemes, \$L=$L\$, $seeds samples"), plot_titlevspan=0.1)

savefig("plots/Graphs/Uniform with Neighbourhood rules over different regiemes.pdf")

println("Saved plot!")