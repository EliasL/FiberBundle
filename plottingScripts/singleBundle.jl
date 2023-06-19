using LaTeXStrings
using Statistics
include("ploting_settings.jl")
include("../burningMan.jl")
include("../support/dataManager.jl")

function breakBundle(b::FB, s::FBS)
    σ = zeros(Float64, b.N)
    for i in 1:b.N
        update_tension!(b)
        find_next_fiber!(b)
        update_storage!(b, s)
        f = b.break_sequence[i]
        tension = b.tension[f]
        σ[i] = tension
        break_fiber!(b)
        resetBundle!(b)
    end
    #average sigma
    every=500
    return σ, map(mean, Iterators.partition(σ, every))
end

function timeBreak(b::FB, s::FBS; real_sigma=false, real_threshold=false,
    force_controlled=false)
    x = [0.0]
    σ = [0.0]
    b.x[argmax(b.x)] = 0.90
    current_x = 0
    last_break_x = 0
    last_break_σ = 0
     # Update bundle stuff
    update_tension!(b)
    find_next_fiber!(b)
    update_storage!(b, s)

    while b.current_step <= b.N
        if real_sigma
            current_σ = sum(current_x*b.σ)/b.N 
        else
            current_σ = current_x * (b.N-b.current_step+1)/b.N
        end
        
        push!(x, current_x)
        if force_controlled
            push!(σ, maximum([current_σ, last_break_σ]))
        else
            push!(σ, current_σ)
        end
        #Check if a fiber will break
        # Curent weakest fiber
        weak = b.break_sequence[b.current_step]
        if real_threshold
            threshold = b.x[weak] / b.σ[weak]
        else
            threshold = b.x[weak]
        end
        if current_x >= threshold 
            push!(x, current_x)
            if force_controlled
                push!(σ,  maximum([current_σ, last_break_σ]))
            else
                push!(σ, current_σ)
            end
            last_break_x = current_x
            if last_break_σ < current_σ
                last_break_σ = current_σ
            end
            current_x = current_σ*0.7 - 0.05*b.current_step/b.N

            # Update bundle stuff
            break_fiber!(b)
            resetBundle!(b)
            update_tension!(b)
            if b.current_step==b.N
                push!(x, last_break_x)
                if force_controlled
                    push!(σ, last_break_σ)
                else
                    push!(σ, 0)
                end
                break 
            else
                find_next_fiber!(b)
            end
            update_storage!(b, s)

        end
        current_x += 1/(1000*b.N)
        
    end
    return x, σ
end
function slowBreak(b::FB, s::FBS; real_sigma=false, real_threshold=false)
    x = [0.0]
    σ = [0.0]
    b.x[argmax(b.x)] = 0.90
    current_x = 0
    last_break_x = 0
    last_break_σ = 0
     # Update bundle stuff
    update_tension!(b)
    find_next_fiber!(b)
    update_storage!(b, s)

    while b.current_step <= b.N
        if real_sigma
            current_σ = sum(current_x*b.σ)/b.N 
        else
            current_σ = current_x * (b.N-b.current_step+1)/b.N
        end
        
        if current_x > last_break_x || current_σ > last_break_σ
            push!(x, current_x)
            push!(σ, current_σ)
        end
        #Check if a fiber will break
        # Curent weakest fiber
        weak = b.break_sequence[b.current_step]
        if real_threshold
            threshold = b.x[weak] / b.σ[weak]
        else
            threshold = b.x[weak]
        end
        if current_x >= threshold 
            push!(x, current_x)
            push!(σ, current_σ)
            last_break_x = current_x
            last_break_σ = current_σ - b.tension[weak]/b.N
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
        current_x += 1/(1000*b.N)
        
    end
    return x, σ
end

function slowBreak2(b::FB, s::FBS; real_sigma=false, real_threshold=false)
    x = [0.0]
    σ = [0.0]
    b.x[argmax(b.x)] = 0.90
    current_x = 0
    last_break_x = 0
    last_break_σ = 0
     # Update bundle stuff
    update_tension!(b)
    find_next_fiber!(b)
    update_storage!(b, s)

    while b.current_step <= b.N
        if real_sigma
            current_σ = sum(current_x*b.σ)/b.N 
        else
            current_σ = current_x * (b.N-b.current_step+1)/b.N
        end
        
        push!(x, current_x)
        push!(σ, current_σ)
        if current_x > last_break_x || current_σ > last_break_σ
        end
        #Check if a fiber will break
        # Curent weakest fiber
        weak = b.break_sequence[b.current_step]
        if real_threshold
            threshold = b.x[weak] / b.σ[weak]
        else
            threshold = b.x[weak]
        end
        if current_x >= threshold 
            push!(x, current_x)
            push!(σ, current_σ)
            last_break_x = current_x
            last_break_σ = current_σ - b.tension[weak]/b.N
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
        current_x += 1/(1000*b.N)
        
    end
    return x, σ
end


function make_plot(b::FB, s::FBS)
    x, σ = slowBreak(b, s, real_sigma=false, real_threshold=false)
    p1 = plot(x, σ, legend=:topleft, title="No Load Sharing ", label="",
            c=:black, xlabel=L"x", ylabel=L"\tilde{σ}",
            ylims=(0, Inf), xlims=(0, 1), size=(300,250))
        x1 = 4100
        x2 = 7115
        plot!([0 0 ; x[x1] x[x2]], [0 0 ; σ[x1] σ[x2]],
        linestyle=:dash, label="", c=:black)
    healBundle!(b)
    x, σ = slowBreak2(b, s, real_sigma=false, real_threshold=true)
    p4 = plot(x, σ, legend=:topleft, title="Equal Load Sharing", label="", size=(300,250), 
            c=:black, xlabel=L"F", ylabel=L"\tilde{σ}", ylims=(0, maximum(σ)*1.0), xlims=(0, Inf))
    healBundle!(b)
    x, σ = slowBreak(b, s, real_sigma=true, real_threshold=true)
    p2 = plot(x, σ, legend=:topleft, title="A", label="", c=:black, 
            xlabel="x", ylabel=L"σ", ylims=(0, maximum(σ)*1.1), xlims=(0, Inf))
    healBundle!(b)
    x, σ = breakBundle(b, s)
    p3 = scatter(1:length(x), x, legend=:topleft, title="B", label="", c=:black, 
            xlabel="k", ylabel=L"σ, x", ylims=(0, maximum(x)*1.1), xlims=(0.5, 9.5),
            xticks=1:9)
    #p3 = make_k_plot(b, s)
#=     healBundle!(b)
    x, σ = timeBreak(b, s, real_sigma=true, real_threshold=true, force_controlled=true)
    plot!((1:length(σ))./length(σ)*10, σ, legend=:topleft, title="A: "* b.nr, label="", c=:black, 
            xlabel="time", ylabel=L"σ, x", ylims=(0, maximum(σ)*1.1), xlims=(0, 11), linestyle=:dash) =#
    savefig(p1, "plots/Graphs/ForceOfSingleBundleWithoutLoadSharing.pdf")
    savefig(p4, "plots/Graphs/ForceOfSingleBundleWithLoadSharing.pdf")
    return p2, p3
end


function make_k_plot(b::FB, s::FBS)
    σ = breakBundle(b, s) 
    x = (1:length(σ)) ./ length(σ)
    p1 = plots(x, σ, legend=:topleft,  label="",
            c=:black, xlabel=L"k/N", ylabel=L"σ",
            alpha=0.1)
    return p1 
end

nr = "ELS"
t = 0.5
L=3
α = 2.0
seed = 5
dist="ConstantAverageUniform"
data_path="newData/"
println("Running...")
b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist)
println(sort(b.x))
println("Making plot...")
ELSPlot = make_plot(b,s)
p=plot(ELSPlot..., size=(300*length(ELSPlot), 300), layout= @layout([ A B;]))
savefig("plots/Graphs/ForceOfSingleBundle.pdf")
#display(p)
#= for nr = ["LLS", "CLS"]
    p = plot()

    for seed= 1:200 
        b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist)
        CLSPlot = make_k_plot(b,s)
        
        #f = load_file(L, α, t, nr, dist, data_path=data_path, average=false)
        #println(maximum(f["most_stressed_fiber/$seed"])*L*L)
        
        p=plot(CLSPlot, size=(500, 300), title=nr, layout= @layout([ A;]))
    end
    display(p)
    savefig("plots/Graphs/$(nr)ForceOfSingleBundle.pdf")
end =#