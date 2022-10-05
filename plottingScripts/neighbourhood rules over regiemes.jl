using Plots
using JLD2
using LaTeXStrings

include("../support/ploting_settings.jl")

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*".jld2"

global_path = "data/"

distributions(t) = [ "t=$t Uniform", "t=$t Uniform SNR", "t=$t Uniform CNR"]
lables = ["something else"]
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
lables = permutedims([d for d in lables])
legend_lables = permutedims(zeros(Int64, length(distributions)))

files = [[file(global_path, L, distribution) for distribution in distributions(t)] for t in (1:9) ./ 10]
seeds = files[1][1]["nr_seeds_used"]
#Files [t] [distribution] [value]

nr_clusters_plot = plot(k_N, [f["average_nr_clusters"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"\#/N", title="Relative number of clusters", legend=false)

largest_cluster_plot = plot(k_N, [f["average_largest_cluster"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"S_{\mathrm{max}}/N", title="Size of largest cluster", legend=false)

largest_perimiter_plot = plot(k_N, [f["average_largest_perimiter"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"H_{\mathrm{max}}/N", title="Length of the longest perimeter", legend=false)

most_stressed_fiber_plot = plot(k_N, [f["average_most_stressed_fiber"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"σ", title="Stress of most stressed fiber", legend=false)

legend_plot = plot(legend_lables, axis=nothing, showaxis = false, grid = false, label=lables, legend=:inside)


l = @layout [
    A B; E; C D
]
plot(nr_clusters_plot, largest_cluster_plot, legend_plot, largest_perimiter_plot, most_stressed_fiber_plot, layout=l,
    plot_title=latexstring("Neighbourhood rules \$t_0=$t, L=$L, $seeds\$ s."), plot_titlevspan=0.1)

savefig("plots/Graphs/Uniform with Neighbourhood rules t₀=$t.pdf")

println("Saved plot!")