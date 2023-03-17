#= using Plots
using Colors
using LaTeXStrings
#using PyPlot
Plots.reset_defaults()
Plots.resetfontsizes()
pyplot()


plot_font = "Comic Sans ms"
default(fontfamily=plot_font, titlefont=font(20, plot_font), foreground_color_legend = nothing, background_color_legend = nothing, 
        grid=false, framestyle=:box, markeralpha=0.0, markerstrokealpha=1, markerstrokewidth=0.7)
Plots.scalefontsizes(1.1)

plot([1,2,3,3], xaxis=L"Test_a", title=L"Test_B", fontfamily="Comic Sans ms")

 =#

using Plots
using LaTeXStrings
pyplot()
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.3)

plot(sort(rand(10)),sort(rand(10)),label="Legend")
plot!(xlabel=L"\textrm{Standard text}(r) / \mathrm{cm^3}")
plot!(ylabel="Same font as everything")
annotate!(0.2,0.8,text("My note",plot_font,12))