using Plots
gr()
Plots.reset_defaults()
plot_font = "Computer Modern"
Plots.resetfontsizes()
default(fontfamily=plot_font, grid=true, legendfontsize=8, markerstrokewidth=0,
        linewidth=2)
Plots.scalefontsizes(1)