using Plots
#gr() # The same issue is also present in the gr backend
pyplot()

function make_plot(X, Y, add_plot=true)

    p = plot()
    if add_plot
        plot!([1,2],[1,2])
        plot!([1,2],[1,2])
        plot!([1,2],[1,2])
    end
    tc = [:red, :black, :orange]
    for (i,s) in enumerate(["text", "with", "color"])
        annotate!(1 + 0.25*i , 1 + 0.5, text(s, tc[i], 18))
    end
    return p
end

plot1 = make_plot([1], [1], true)
plot2 = make_plot([1], [1], false)
plot(plot1, plot2)

