using Plots
using Interpolations

include("ploting_settings.jl")

y=[
     8  8 8;
    10 10 10;
    10 2 6;
    12 2 7;
    2  2 2;
    2  2 2;
]./2
x = repeat([0, 1, 1, 2, 2, 3].*2, 1, 2)
p = plot(x, y, size=(250,200), style=[:solid :solid :dash], label=[L"A" L"B" ""],
    xlabel=L"k", ylabel=L"\sigma")
savefig(p, "plots/Graphs/illustrativeExample.pdf")