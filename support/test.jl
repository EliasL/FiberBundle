test1(x) = sum(x[x .!= Inf])

function test2(x)
    s = 0
    for i in eachindex(x)
        @inbounds s += ifelse(isinf(x[i]), 0, x[i])
    end
    return 0
end

x = fill(Inf, 10^8)
x = vcat(x, [1,2,3])

@time test1(x)
@time test1(x)
@time test2(x)
@time test2(x)
