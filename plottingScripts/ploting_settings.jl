using Plots
using Colors
Plots.reset_defaults()
Plots.resetfontsizes()
gr()
plot_font = "Computer Modern"
default(fontfamily=plot_font, grid=true, legendfontsize=8, markerstrokewidth=0, markersize=2.4,
        linewidth=2, foreground_color_legend = nothing, background_color_legend = RGBA(1.0,1.0,1.0,0.3))
Plots.scalefontsizes(1)