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
    return m
end

function generate_neighbours()
    nr_neighbours = 256
    neighbours = [zeros(Int, (3, 3)) for i in 1:nr_neighbours]
    for i in 1:nr_neighbours
        n = reverse(parse.(Int, split(string(i-1, base=2, pad=9), "")))
        neighbours[i] = reshape_neighbours(n)
    end
    return neighbours
end


function remove_symmetries(n)
    n = copy(n)
    for grid in n
        syms = [grid, rotr90(grid), rot180(grid), rotl90(grid), reverse(grid, dims=1), reverse(grid, dims=2), rotr90(reverse(grid, dims=1)), rotr90(reverse(grid, dims=2))]
        duplicates = findall(x->x in syms,n)[2:end]
        deleteat!(n, duplicates)
    end
    return n
end

function compute_strength(m::Matrix)
    # Here are the strength values
    # 121
    # 2_2
    # 121
    strength = abs.(sum([m[2,1], m[1,2], m[3,2], m[2,3]]))*2
    strength += abs.(sum([m[1,1], m[3,1], m[1,3], m[3,3]]))

    # Rotation symetry
    if m == rot180(m)
        strength += 0.5
    end
    return strength
end

function neighbourhoodToInt(m::AbstractArray{Int64})
    id = 1# Not zero because of zero indexing in julia
    for (i, alive) in enumerate(m)
        if alive<0
            if i > 5
                id +=2^(i-2)
            elseif i < 5
                id +=2^(i-1)
            end
        end
    end
    return id
end

neighbourhoodStrengths = zeros(Float64, 256)
for m in generate_neighbours().*-1
    i = neighbourhoodToInt(m)
    s = compute_strength(m)
    neighbourhoodStrengths[i] = s
end