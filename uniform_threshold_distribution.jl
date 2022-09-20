using Plots
using JLD2
using LaTeXStrings

include("ploting_settings.jl")

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*".jld2"

global_path = "data/"
distribution = "Uniform"

function file()
    return load(full_name(global_path, L, distribution))
end

L = [32, 64, 128, 256]
N = L.*L
k_N = [1:n for n in N]./N
lables = permutedims([latexstring("\$L={$(l)}\$") for l in L])

files = [load(full_name(global_path, L, distribution)) for L in [32, 64, 128, 256]]
seeds = files[1]["seeds_used"]

nr_clusters_plot = plot(k_N, [f["average_nr_clusters"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"\#/N", title="Relative number of clusters")

largest_cluster_plot = plot(k_N, [f["average_largest_cluster"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"S_{\mathrm{max}}/N", title="Size of largest cluster", legend=:topleft)

largest_perimiter_plot = plot(k_N, [f["average_largest_perimiter"] for f in files], label = lables,
                    xlabel=L"k/N", ylabel=L"H_{\mathrm{max}}/N", title="Length of the longest perimeter")

plot(nr_clusters_plot, largest_cluster_plot, largest_perimiter_plot, layout=(2,2),
    plot_title="Uniform distribution, $seeds samples", plot_titlevspan=0.1)
savefig("plots/Uniform.pdf")
println("Saved plot!")