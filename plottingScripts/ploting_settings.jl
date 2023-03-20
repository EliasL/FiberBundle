using Plots
using Colors
#using PyPlot
Plots.reset_defaults()
Plots.resetfontsizes()
pyplot()


plot_font = "Computer Modern"
default(#= fontfamily=plot_font, =#foreground_color_legend = nothing,
        background_color_legend = nothing, grid=false, framestyle=:box,
        markeralpha=0.0, markerstrokealpha=1, markerstrokewidth=0.7)
Plots.scalefontsizes(1.1)