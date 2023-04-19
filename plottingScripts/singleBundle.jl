using LaTeXStrings
include("ploting_settings.jl")
include("../burningMan.jl")
include("../support/dataManager.jl")

function breakBundle(b::FB, s::FBS)
    x = zeros(Float64, 2*b.N-1)
    σ = zeros(Float64, 2*b.N-1)
    for i in 1:b.N-1
        update_tension!(b)
        find_next_fiber!(b)
        update_storage!(b, s)
        f = b.break_sequence[i]
        X = b.x[f]#/b.σ[f] 
        x[2*i] = X
        x[2*i+1] = x[2*i]
        σ[2*i] = b.max_σ
        σ[2*i-1] = σ[2*i] - b.max_σ/X * (X-x[2*i-1])
        break_fiber!(b)
        
        resetBundle!(b)
    end
    return x, σ
end

function slowBreak(b::FB, s::FBS)
    x = [0.0]
    σ = [0.0]
    current_x = 0
    last_break_x = 0
     # Update bundle stuff
    update_tension!(b)
    find_next_fiber!(b)
    update_storage!(b, s)

    while b.current_step <= b.N
        current_σ = (b.N-b.current_step)*current_x
        if current_x >= last_break_x
            push!(x, current_x)
            push!(σ, current_σ)
        end
        #Check if a fiber will break
        # Curent weakest fiber
        weak = b.break_sequence[b.current_step]
        threshold = b.x[weak]
        real_threshold = threshold / b.σ[weak]
        if current_x >= real_threshold 
            # Update bundle stuff
            break_fiber!(b)
            resetBundle!(b)
            update_tension!(b)
            if b.current_step==b.N
                push!(x, 1)
                push!(σ, 0)
                break 
            else
                find_next_fiber!(b)
            end
            update_storage!(b, s)

            last_break_x = current_x
            current_x = 0
        end
        current_x += 1/(10*b.N)
        
    end
    return x, σ
end

function make_plot(b::FB, s::FBS)
    x, σ = breakBundle(b, s)
    return plot(x, σ, legend=:topleft, title="Not " * b.nr, label="", c=:black, xlabel="sort of "*L"x", ylabel="almost " * L"σ")
end


nr = "ELS"
t = 0.5
L=4
α = 2.0
seed = 1
dist="ConstantAverageUniform"
data_path="newData/"

b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist)
ELSPlot = make_plot(b,s)

nr = "LLS"
b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist)
LLSPlot = make_plot(b,s)

nr = "CLS"
b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist)
CLSPlot = make_plot(b,s)

#f = load_file(L, α, t, nr, dist, data_path=data_path, average=false)
#println(maximum(f["most_stressed_fiber/$seed"])*L*L)

plot(ELSPlot, size=(300, 300), layout= @layout([ A;]), link=:y)
savefig("plots/Graphs/ForceOfSingleBundle.pdf")