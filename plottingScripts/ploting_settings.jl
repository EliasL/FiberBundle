using Plots
using Colors
Plots.reset_defaults()
Plots.resetfontsizes()
gr()
plot_font = "Computer Modern"
default(fontfamily=plot_font, foreground_color_legend = nothing, background_color_legend = RGBA(1.0,1.0,1.0,0.3),
        grid=false)
Plots.scalefontsizes(1.3)