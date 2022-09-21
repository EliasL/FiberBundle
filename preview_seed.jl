using Measures
using JLD2
using LaTeXStrings

using CairoMakie
using FileIO 

#include("ploting_settings.jl")

full_name(global_path, L, distribution) = global_path*distribution*"/"*distribution*string(L)*".jld2"
load_file() = load(full_name(global_path, L, distribution))

global_path = "data/"
distribution = "Uniform"
L = 64
seed = 1
file = load_file()


f = Figure(backgroundcolor = RGBAf0(0.98, 0.98, 0.98))
ga = f[1, 1] = GridLayout()
axies = [Axis(ga[i,j]) for i=1:3, j=1:3]

for (i, axis) in enumerate(axies)
    g = reshape(file["sample_states/$seed"][i, :], (L, L))
    fig, ax, pltobj = heatmap(g)
end

f
# grids = [heatmap(reshape(f["sample_states/$seed"][i, :], (L, L)), colorbar=false, aspect_ratio=:equal, axis=([], false), clim=(0, 1), title=L"k/N = 0."*"$i") for i in 1:9]

#  plot(grids..., layout=(3,3),
#     plot_title="$distribution distribution, "*L"L="*"$L", plot_titlevspan=0.1)
# savefig("plots/Uniform.pdf")
# println("Saved plot!")