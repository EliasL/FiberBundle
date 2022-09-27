using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*".jld2"

global_path = "data/"
distributions = ["Uniform with Neighbourhood rules", "Uniform"]
desired_data = ["average_nr_clusters", "average_largest_cluster", "average_largest_perimiter", "nr_seeds_used"]

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

L = 256
N = L.*L
k_N = [1:n for n in N]./N
lables = permutedims([d for d in distributions])

files = [file(global_path, L, distribution) for distribution in distributions]
seeds = files[1]["nr_seeds_used"]

nr_clusters_plot = plot(k_N, [f["average_nr_clusters"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"\#/N", title="Relative number of clusters", legend=false)

largest_cluster_plot = plot(k_N, [f["average_largest_cluster"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"S_{\mathrm{max}}/N", title="Size of largest cluster", legend=false)

largest_perimiter_plot = plot(k_N, [f["average_largest_perimiter"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"H_{\mathrm{max}}/N", title="Length of the longest perimeter", legend=:outerright)


l = @layout [
    a{0.5w} b{0.5w}
    c{0.5h} 
]
plot(nr_clusters_plot, largest_cluster_plot, largest_perimiter_plot, layout=l,
    plot_title="Neighbourhood rules, $seeds samples, L=$L", plot_titlevspan=0.1)
savefig("plots/Uniform with Neighbourhood rules.pdf")
println("Saved plot!")