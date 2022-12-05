using Plots
gr()
default(markershape=:circle)
a = plot([1 1; 2 2], [1 2; 2 3], linestyle=[:solid :dash],
    labels=["solid" "dashed"], size=(2000,1000))
