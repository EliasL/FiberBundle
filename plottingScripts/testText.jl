using Plots

pyplot()
Plots.reset_defaults()
default(framestyle=:box)

p = plot(100rand(10), c=:green, framestyle=:box)
plot!(p, rand(10), c=:blue, inset=bbox(0, 0, 0.4, 0.4, :center), subplot=2)
plot!(twinx(p[2]), 10rand(10), c=:red)
plot!(p[2], framestyle=:box)