using Plots
using Colors
#using PyPlot
Plots.reset_defaults()
Plots.resetfontsizes()
pyplot()
plot_font = "Comic Sans MS"#"Computer Modern"
default(fontfamily=plot_font, foreground_color_legend = nothing, background_color_legend = nothing,
        grid=false, framestyle=:box,)
Plots.scalefontsizes(1.1)