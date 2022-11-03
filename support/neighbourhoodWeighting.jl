using BenchmarkTools

function reshape_neighbours(l)
    m = reshape(l, (3,3))
    # Shifts the numbers to make room for
    # the center piece
    m[3,3] = m[2,3]
    m[2,3] = m[1,3]
    m[1,3] = m[3,2]
    m[3,2] = m[2,2]
    m[2,2] = 2
    return vec(m)
end

function generate_neighbours(;return_matrix = false)
    if return_matrix
        l = 8
    else
        l = 9
    end
    nr_neighbours = 256
    neighbours = [zeros(Int, l) for i in 1:nr_neighbours]
    for i in 1:nr_neighbours
        neighbours[i] = reverse(parse.(Int, split(string(i-1, base=2, pad=9), "")))
        if return_matrix
            neighbours[i] = reshape_neighbours(neighbours[i])
        end
    end
    return neighbours
end


function remove_symmetries(n)
    n = [reshape(v, (3,3)) for v in n]
    for grid in n
        syms = [grid, rotr90(grid), rot180(grid), rotl90(grid), reverse(grid, dims=1), reverse(grid, dims=2), rotr90(reverse(grid, dims=1)), rotr90(reverse(grid, dims=2))]
        duplicates = findall(x->x in syms,n)[2:end]
        deleteat!(n, duplicates)
    end
    return vec.(n)
end

function compute_strength(m::Vector{Int64})
    m = reshape(m, (3,3))
    # Here are the strength values
    # 121
    # 2_2
    # 121
    strength = abs.(sum([m[2,1], m[1,2], m[3,2], m[2,3]]))*4
    strength += abs.(sum([m[1,1], m[3,1], m[1,3], m[3,3]]))*2

    # Rotation symetry
    if m == rot180(m)
        strength += 1
    end
    return strength
end

#https://discourse.julialang.org/t/parse-an-array-of-bits-bitarray-to-an-integer/42361/24
function arr_to_int(a, val = 0)
    v = 128
    # Check for negative numbers. Negative status means alive
    for i in [a[1]<0, a[2]<0, a[3]<0,
              a[4]<0,         a[5]<0,
              a[6]<0, a[7]<0, a[8]<0]
        val += v*i
        v >>= 1
    end
    return val
end

function neighbourhoodToStrength(a::AbstractVector{Int64})
    # Add one because of 1 indexing
    i =  arr_to_int(a, 1)
    return neighbourhoodStrengths[i]
end

neighbourhoodStrengths = zeros(Int64, 256)
for m in generate_neighbours().*-1
    println(m)
    println(neighbourhoodToStrength(m))
    i = neighbourhoodToStrength(m)
    str = compute_strength(m)
    @assert neighbourhoodStrengths[i] == 0 "This value is already set! arr_to_int does not work!"
    neighbourhoodStrengths[i] = str
end

