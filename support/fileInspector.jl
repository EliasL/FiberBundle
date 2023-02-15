

include("dataManager.jl")
include("../plottingScripts/showBundle.jl")

function show_bundle(l, t, nr, α, seed, progression=0.0)
    settings = make_settings(l, t, nr, α)
    b = get_bundles_from_settings(settings, seeds=seed, spanning=false, progression=progression)
    p = plot_fb(b)
end


function print_d(d)
    if d isa AbstractArray
        s = size(d)
        if length(s) == 1
            print("l: $(s[1]), ")
        else
            print("s: $(size(d)), ")
        end
        println(round(first(d), sigdigits=2), " ... ", round(last(d), sigdigits=2))
    elseif d isa Number
        println(round(d, digits=3))
    else
        println(typeof(d))
    end
end

L=32
a=2.0
nr="CLS"
t=0.9
seed = 9999
f_path = get_file_path(L, a, t, nr, average=false)
f = load(f_path)
println("Loaded $f_path")

if haskey(f, "average_nr_clusters")
    println("Averages: ")
    for key in averaged_data_keys
        print("$key => ")
        d = f["average_$key"]
        print_d(d)
    end
else
    println("Seed $seed: ")
    for key in data_keys
        #display(f)
        print("$key => ")
        d = f["$key/$seed"]
        print_d(d)
    end
end
seeds = f["seeds_used"]
except = setdiff(first(seeds):last(seeds), seeds)
println("Seeds: $(first(seeds)) - $(last(seeds))" * (isempty(except) ? "" : " except $(join(seeds, " , "))"))
show_bundle(L, t, nr, α, seed, 0.875)