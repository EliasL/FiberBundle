using Plots
using LaTeXStrings
gr()
nr="LLS"
slope_h = 2
h_err = 3
plot(title=latexstring("$nr: \$D_h=$slope_h\$")*" ± "*"$h_err")
