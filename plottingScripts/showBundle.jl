using Plots

include("../burningMan.jl")

function show_fb(b::FB)
    show_array(b.status, b.L)
end

function show_array(a::AbstractArray, L=nothing)
    if L===nothing
        L = round(Int, sqrt(length(a)))
    end
    m = reshape(a, (L, L))
    show_matrix(m)
end
function show_matrix(m::AbstractMatrix)
    h = heatmap(m, c=:glasbey_category10_n256, legend=:none, showaxis = false, ticks=false)
    plot(h)
end
