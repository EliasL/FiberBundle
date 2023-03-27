using Random


function get_fixed_average_uniform_distribution(t0)
    return x(n) = 0.5 .+ (-t0 .+ 2*rand(Float64, n)*t0)
end


function get_uniform_distribution(t0)
    return x(n) = 1 .- (rand(Float64, n) .* (1-t0))
end

function uniform(n::Int64)
    return rand(Float64, n)
end

