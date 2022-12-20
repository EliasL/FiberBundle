using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*".jld2"

global_path = "data/"
distribution = "t=0.0 Uniform"

function file()
    return load(full_name(global_path, L, distribution))
end

L = [32, 64, 128]
N = L.*L
k_N = [1:n for n in N]./N
lables = permutedims([latexstring("\$L={$(l)}\$") for l in L])
legend_lables = permutedims(zeros(Int64, length(L)))

files = [load(full_name(global_path, l, distribution)) for l in L]
seeds = files[1]["nr_seeds_used"]

nr_clusters_plot = plot(k_N, [f["average_nr_clusters"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"\#/N", title="Relative number of clusters")

largest_cluster_plot = plot(k_N, [f["average_largest_cluster"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"S_{\mathrm{max}}/N", title="Size of largest cluster", legend=:topleft)

largest_perimiter_plot = plot(k_N, [f["average_largest_perimiter"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"H_{\mathrm{max}}/N", title="Length of the longest perimeter")

most_stressed_fiber_plot = plot(k_N, [f["average_most_stressed_fiber"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"σ", title="Stress of most stressed fiber", legend=false)

legend_plot = plot(legend_lables, axis=nothing, showaxis = false, grid = false, label=lables, legend=:inside)


l = @layout [
    A B; E; C D
]
plot(nr_clusters_plot, largest_cluster_plot, legend_plot, largest_perimiter_plot, most_stressed_fiber_plot, layout=l,
    plot_title="Uniform distribution, $seeds samples", plot_titlevspan=0.1)
    
savefig("plots/Graphs/Uniform.pdf")
println("Saved plot!")