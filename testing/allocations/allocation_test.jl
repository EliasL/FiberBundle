using StaticArrays


Base.@kwdef mutable struct FB{F<:AbstractFloat, I<:Integer}
    L::I
    N::I = L*L
    x::Vector{F} =ones(F, N)
    σ::Vector{F} = ones(F, N) # Relative tension
    tension::Vector{F} = zeros(F, N) 
    max_σ::F = 0.0
    break_sequence::Vector{I} = zeros(I, N) 
end


function findNextFiber!(b::FB)
        b.tension[1] = b.σ[1] / b.x[1]
        b.break_sequence[1] = argmax(b.tension)
        b.max_σ = b.tension[b.break_sequence[1]]
end

function dict_make_FB(L)
    settings = Dict(
        "L" => L,
    )
    L = Int64(settings["L"])
    return FB{Float64, Int64}(L=L)
end

function make_FB(L)
    return FB{Float64, Int64}(L=L)
end

function test()
    L=8
    b1 = make_FB(L)
    b2 = dict_make_FB(L)

    findNextFiber!(b1) #Compile

    println("Does not have allocation")
    for _ in 1:2
        @time findNextFiber!(b1)
    end
    println("Has allocation")
    for _ in 1:2
        @time findNextFiber!(b2)
    end

end
test()