using LaTeXStrings
include("ploting_settings.jl")
include("../burningMan.jl")
include("../support/dataManager.jl")

function breakBundle(b::FB, s::FBS, real_values=false)
    x = zeros(Float64, 2*b.N+1)
    σ = zeros(Float64, 2*b.N+1)
    b.x[argmax(b.x)] = 0.92
    for i in 1:b.N
        update_tension!(b)
        find_next_fiber!(b)
        update_storage!(b, s)
        f = b.break_sequence[i]
        if real_values #Not quite working
            X = b.x[f]/b.σ[f] 
            tension = b.tension[f]*(b.N-b.current_step+1)
            x[2*i] = X
            x[2*i+1] = x[2*i]
            σ[2*i] = tension
            σ[2*i-1] = tension - b.max_σ
        else
            X = b.x[f] 
            tension = X * (b.N-b.current_step+1)
            x[2*i] = X
            x[2*i+1] = x[2*i]
            σ[2*i] = tension
            σ[2*i+1] = tension-X
        end
        break_fiber!(b)
        
        resetBundle!(b)
    end
    return x, σ
end

function slowBreak(b::FB, s::FBS, real_sigma=false)
    x = [0.0]
    σ = [0.0]
    b.x[argmax(b.x)] = 0.92
    current_x = 0
    last_break_x = 0
    last_break_σ = 0
     # Update bundle stuff
    update_tension!(b)
    find_next_fiber!(b)
    update_storage!(b, s)

    while b.current_step <= b.N
        if real_sigma
            current_σ = sum(current_x*b.σ) 
        else
            current_σ = current_x * (b.N-b.current_step+1)
        end
        
        if current_x >= last_break_x || current_σ >= last_break_σ
            push!(x, current_x)
            push!(σ, current_σ)
        end
        #Check if a fiber will break
        # Curent weakest fiber
        weak = b.break_sequence[b.current_step]
        threshold = b.x[weak]
        real_threshold = threshold / b.σ[weak]
        if current_x >= real_threshold 
            push!(x, current_x)
            push!(σ, current_σ)
            last_break_x = current_x
            last_break_σ = current_σ - b.tension[weak]*1.8
            current_x = 0

            # Update bundle stuff
            break_fiber!(b)
            resetBundle!(b)
            update_tension!(b)
            if b.current_step==b.N
                push!(x, last_break_x)
                push!(σ, 0)
                break 
            else
                find_next_fiber!(b)
            end
            update_storage!(b, s)

        end
        current_x += 1/(50*b.N)
        
    end
    return x, σ
end

function make_plot(b::FB, s::FBS)
    x, σ = breakBundle(b, s)
    p1 = plot(x, σ, legend=:topleft, title="No Load Sharing ", label="", c=:black, xlabel="sort of "*L"x", ylabel="slightly " * L"σ")
    healBundle!(b)
    x, σ = slowBreak(b, s)
    p2 = plot(x, σ, legend=:topleft, title="Almost " * b.nr, label="", c=:black, xlabel=L"x", ylabel=L"\tilde{σ}")
    healBundle!(b)
    x, σ = slowBreak(b, s, true)
    p3 = plot(x, σ, legend=:topleft, title=b.nr, label="", c=:black, xlabel=L"x", ylabel=L"σ")
    return p3, p2, p1
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

plot(ELSPlot..., size=(900, 300), layout= @layout([ A B C;]))
savefig("plots/Graphs/ForceOfSingleBundle.pdf")