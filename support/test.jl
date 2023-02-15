using Plots

# Generate some data
x = 1:10
y1 = rand(10)
y2 = rand(10)

# Create the first plot with left y-axis
p1 = plot(x, y1, label="Plot 1", xlabel="", ylabel="Y-axis", xticks=x)

# Create the second plot with right y-axis and no label or ticks
p2 = plot(x, y2, label="Plot 2", xlabel="", xticks=[], showaxis=[false, true], tickfont=font(2))

plot!(p2, p1, layout = (2, 1), link=:x)
