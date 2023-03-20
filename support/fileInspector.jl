

include("dataManager.jl")
include("../plottingScripts/showBundle.jl")

function show_bundle(l, t, nr, α, seed, dist; progression=0.0, critical=false)
    settings = make_settings(l, t, nr, α, dist=dist)
    b = get_bundles_from_settings(settings, seeds=seed, spanning=false, critical=critical, progression=progression)
    p = plot_fb(b, cm_shift=true)
end


function print_d(d)
    if d isa AbstractArray
        s = size(d)
        if length(s) == 1
            println("\tlength: $(s[1]), ")
        else
            println("\tsize: $(size(d)), ")
        end
        println("\targmax: $(argmax(d))")
        println("\tmax: $(maximum(d))")
        println("\taverage: $(sum(d)/length(d))")
        println("\t ",round(first(d), sigdigits=2), " ... ", round(last(d), sigdigits=2))
    elseif d isa Number
        println("\t", round(d, digits=3))
    else
        println("\t", typeof(d))
    end
end

function print_file(f_path, seed=seed)
    jldopen(f_path, "r") do f
        println("Loaded $f_path")
    
        if haskey(f, "average_nr_clusters")
            println("Averages: ")
            for key in averaged_data_keys
                println("$key => ")
                d = f["average_$key"]
                print_d(d)
            end
        else
            println("Seed $seed: ")
            for key in data_keys
                #display(f)
                println("$key => ")
                d = f["$key/$seed"]
                print_d(d)
            end
        end
        seeds = f["seeds_used"]
        except = setdiff(first(seeds):last(seeds), seeds)
        println("Seeds: $(first(seeds)) - $(last(seeds))" * (isempty(except) ? "" : " except $(join(seeds, " , "))"))
        #println(f["largest_cluster/995"])
    end
end

function investigate(f_path, seed)
    jldopen(f_path, "r") do f
        for seed in [8]
            d = f["most_stressed_fiber/$seed"]
            if argmax(d) != 1
                print_file(f_path, seed)
                show_bundle(L, t, nr, α, seed, critical=true)
                settings = make_settings(L, t, nr, α, dist)
                b1 = get_bundles_from_settings(settings, seeds=seed, spanning=false, critical=false, progression=1/(L*L))
                b2 = get_bundles_from_settings(settings, seeds=seed, spanning=false, critical=false, progression=2/(L*L))
                b3 = get_bundles_from_settings(settings, seeds=seed, spanning=false, critical=false, progression=3/(L*L))
                plot_fb(b1)
                println(b1.σ)
                println(b1.x)
                plot_fb(b2)
                plot_fb(b3)
                
    
            end
        end
    end
end

L=16
a=2.0
nr="CLS"
t=0.3
dist = "ConstantAverageUniform"
f_path = get_file_path(L, a, t, nr, dist, average=false)

investigate2(f_path,0)
#print_file(f_path)

#= for seed in 995:999
    show_bundle(L, t, nr, α, seed, critical=true)
end =#