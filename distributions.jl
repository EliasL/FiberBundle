using Random

function get_uniform_distribution(t₀)
    return x(n) = 1 .- (rand(Float64, n) .* (1-t₀))
end

function uniform(n::Int64)
    return rand(Float64, n)
end