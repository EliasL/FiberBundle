include("ploting_settings.jl")
include("../burningMan.jl")
include("../support/dataManager.jl")

function forceOnBundle(b::FB, s::FBS)
    x = zeros(Float64, b.N)
    σ = zeros(Float64, b.N)
    for i in 1:b.N
        update_tension!(b)
        find_next_fiber!(b)
        update_storage!(b, s)
        f = b.break_sequence[i]
        X = b.x[f]#/b.σ[f] 
        x[i] = i
        #x[2*i+1] = x[2*i]
        σ[i] = (b.N-i) *b.max_σ
        #σ[2*i+1] = (b.N-i-1) * b.max_σ
        break_fiber!(b)
        
        resetBundle!(b)
    end
    p = plot(x,σ, legend=:topleft, title=b.nr)
    return p
end

function slowBreak(b::FB, s::FBS)
    x = []
    σ = []
    current_x = 0
     # Update bundle stuff
    update_tension!(b)
    find_next_fiber!(b)
    update_storage!(b, s)

    while b.current_step < b.N
        current_σ = (b.N-b.current_step)*current_x
        push!(x, current_x)
        push!(σ, current_σ)
        #Check if a fiber will break
        # Curent weakest fiber
        weak = b.break_sequence[b.current_step]
        threshold = b.x[weak]
        real_threshold = threshold / b.σ[weak]
        if current_x > threshold 
            # Update bundle stuff
            break_fiber!(b)
            resetBundle!(b)
            update_tension!(b)
            find_next_fiber!(b)
            update_storage!(b, s)
            current_x = 0
        end
        current_x += 1/(2*b.N)
        
    end
    p = plot(x,σ, legend=:topleft, title=b.nr)
    return p
end


nr = "ELS"
t = 0.5
L=6
α = 2.0
seed = 1
dist="ConstantAverageUniform"
data_path="newData/"

b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist)
ELSPlot = forceOnBundle(b,s)

nr = "LLS"
b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist)
LLSPlot = forceOnBundle(b,s)
plot(ELSPlot, LLSPlot)
savefig("plots/Graphs/ForceOfSingleBundle.pdf")