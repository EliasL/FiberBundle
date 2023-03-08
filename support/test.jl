using Plots

# Generatesome data
x = ([1:4 , 2:2:8]./8)./[3,9].^0.6333
y2 = ([[0.1, 0.1, 0.1, 0.2], [0.2, 0.2, 0.2, 0.4]]./0.4)./[3,9].^0.63

# Create the second plot with right y-axis and no label or ticks
plot(x, y2, label="Plot 2", xlabel="", yaxis=:log, xaxis=:log)
