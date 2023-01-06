using Plots
using Colors
gr()

Plots.reset_defaults()
plot_font = "Computer Modern"
Plots.resetfontsizes()
default(fontfamily=plot_font, grid=true, legendfontsize=8, markerstrokewidth=0,
        linewidth=2, foreground_color_legend = nothing, background_color_legend = RGBA(1.0,1.0,1.0,0.3))
Plots.scalefontsizes(1)