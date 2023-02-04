using JLD2
using CodecLz4
using Logging
using ProgressMeter
using Statistics

include("logingLevels.jl")
include("../burningMan.jl")

averaged_data_keys = [
                    "simulation_time",
                    "nr_clusters", 
                    "largest_cluster",
                    "largest_perimiter",
                    "most_stressed_fiber",
                    "spanning_cluster_size",
                    "spanning_cluster_perimiter",
                    "spanning_cluster_step",
                    ]
seed_specific_keys = [
                    "last_step",
                    "spanning_cluster_step",
                    "break_sequence",
                    ]
data_keys = Set(vcat(averaged_data_keys, seed_specific_keys))


function make_settings(L::Int64, t::Float64, nr::String, α::Float64, path::String="data/", dist::String="Uniform")
    if nr=="LLS" || nr=="ELS"
        α = 0.0
    end
    settings = Dict(
        "dist" => dist,
        "L" => L,
        "t" => t,
        "nr" => nr,
        "a" => α,
    )
    add_name_and_path_to_setting!(path, dist, settings)

    if !isdir(settings["path"])
        mkpath(settings["path"])
    end
    return settings
end

function add_name_and_path_to_setting!(path, dist, settings) 
    settings["name"] = get_setting_name(settings)
    settings["path"] = path*dist*"/"*settings["name"]*"/"
end

function display_setting(s)
    return "$(s["dist"]) $(s["nr"]) L:$(s["L"]) α:$(s["a"]) t:$(s["t"])"    
end

function get_setting_name(settings)
    name = ""
    excluded = ["name", "path"]
    # We want to sort so that the naming is consistant
    for (setting, value) in sort(collect(settings), by=x->lowercase(x[1]))
        if setting in excluded
            continue
        else
            name *= "$setting=$value "
        end
    end
    return chop(name)
end

function get_file_name(L, α, t, NR, dist="Uniform"; data_path="data/", seed=-1, average=true)
    s = make_settings(L, t, NR, α, data_path, dist)
    return get_file_name(s, seed, average)
end

function get_file_name(settings, seed::Int=-1, average=false)
    @assert seed>=-1 "We don't want to use negative seeds since -1 is special here"
    s = settings
    if average
        return "$(s["path"])$(s["name"]).jld2"
    elseif seed == -1
        return "$(s["path"])$(s["name"])_bulk.jld2"
    else
        return "$(s["path"])$(s["name"]) s=$(seed)_bulk.jld2"
    end
end

function make_get_name(settings)
    return f(seed=-1; average=false) = get_file_name(settings, seed, average)
end

function expand_file(settings, overwritten_seeds::AbstractArray=Vector{Int64}([]), force_overwrite=false)
    get_name_fun = make_get_name(settings)
    condensed_file_name = get_name_fun()
    jldopen(condensed_file_name, "r") do file

        seeds = file["seeds_used"]
        for seed::Int64 in seeds
            if seed in overwritten_seeds
                # If the seed is going to be overwritten, we don't open this file
                # and let it be generated from scratch
                continue
            end
            # Make sure we are not overwriting already existing files
            if !isfile(get_name_fun(seed)) || force_overwrite
                # Write out the file
                jldopen(get_name_fun(seed), "w") do s_file
                    for key in data_keys
                        if haskey(file, "$key/$seed")
                            s_file[key] = file["$key/$seed"]
                        end
                    end
                end # Close seed file
            end
        end
    end # Close compact file
end

function get_seeds_in_file(settings, average=false)
    file_name = get_file_name(settings, -1, average)
    if isfile(file_name)
        jldopen(file_name, "r") do file
            if haskey(file, "seeds_used")
                return file["seeds_used"]
            else
                return []
            end
        end
    else
        return []
    end
end

function get_min_steps_in_files(settings)
    file_name = get_file_name(settings, -1, false)
    if isfile(file_name)
        jldopen(file_name, "r") do file
            if haskey(file, "seeds_used")
                seeds = file["seeds_used"]
                return minimum([file["last_step/$s"] for s in seeds])
            else
                return 0
            end
        end
    else
        return 0
    end
end

function condense_files(settings, requested_seeds::AbstractArray; remove_files=true, keep_existing_seeds=true)
    get_name_fun = make_get_name(settings)
    condensed_file_name = get_name_fun()
    averaged_file_name = get_name_fun(average=true)
    existing_seeds = get_seeds_in_file(settings, )
    if length(existing_seeds) > 0
        expand_file(settings)
    end
    
    seeds = requested_seeds
    if keep_existing_seeds
        seeds = union(existing_seeds, requested_seeds)
    end
    unused_seeds = setdiff(existing_seeds, requested_seeds)

    nr_seeds = length(seeds)

    averages = Dict()

    if settings["L"] > 512
        compress = LZ4FrameCompressor()
    else
        compress = false
    end
    compress = false

    jldopen(condensed_file_name*".temp", "w", compress = compress) do condensed_file
    jldopen(averaged_file_name*".temp", "w", compress = compress) do averaged_file
        averaged_file["seeds_used"] = seeds
        averaged_file["nr_seeds_used"] = length(seeds)
        condensed_file["seeds_used"] = seeds
        condensed_file["nr_seeds_used"] = length(seeds)
        for seed in seeds
            seed_file_name = get_name_fun(seed)
            jldopen(seed_file_name, "r+") do s_file
                for key in data_keys
                    if !haskey(s_file, key)
                        value = 0
                        @warn "$key not found in $(seed_file_name)!"
                    else
                        value = s_file[key]
                    end
                    if key in averaged_data_keys
                        # If the key doesn't exist, make it
                        if !haskey(averages, key)
                                averages[key] = []
                        end
                        push!(averages[key], value)
                    end
                    condensed_file["$key/$seed"] = value
                end
            end
        end
        for key in averaged_data_keys
            m = mean(averages[key])
            s = std(averages[key], mean=m)
            averaged_file["average_$key"] = m
            averaged_file["std_$key"] = s
        end
    end # Averaged file
    end # Condensed file
    if remove_files
        for seed in seeds
            rm(get_name_fun(seed))
        end
    end
    if !keep_existing_seeds
        for seed in unused_seeds
            rm(get_name_fun(seed))
        end
    end
    mv(condensed_file_name*".temp", condensed_file_name, force=true)
    mv(averaged_file_name*".temp", averaged_file_name, force=true)
end

function get_missing_seeds(settings, requested_seeds)

    get_name_fun = make_get_name(settings)
    condensed_file_name = get_name_fun()
    if isfile(condensed_file_name)
        jldopen(condensed_file_name, "r") do existing_data
            existing_seeds = existing_data["seeds_used"]
            return Vector{Int64}(setdiff(requested_seeds, existing_seeds))
        end
    else
        return requested_seeds
    end
end

function prepare_run(settings, requested_seeds::AbstractArray, overwrite=false)
    search_for_loose_files(settings)

    get_name_fun = make_get_name(settings)
    condensed_file_name = get_name_fun()
    # Get missing seeds
    if overwrite
        missing_seeds = requested_seeds
    else
        missing_seeds = get_missing_seeds(settings, requested_seeds)
    end

    # Print some info
    nr_found_files = length(requested_seeds) - length(missing_seeds)
    if nr_found_files>0 && nr_found_files!=length(requested_seeds)
        @logmsg settingLog "$nr_found_files files were found and skipped."
    end

    # Here we expand the existing datafile 
    if isfile(condensed_file_name)
        if length(missing_seeds) == 0
            @logmsg settingLog "All seeds already exist"
        else
            expand_file(settings, missing_seeds)
        end
    end

    return missing_seeds
end

function clean_after_run(settings, requested_seeds::AbstractArray)
    condense_files(settings, requested_seeds, remove_files=true)
end

function search_for_loose_files(settings)
    
    @logmsg settingLog "Searching for loose files... "
    files = readdir(settings["path"])
    # Find distribution with lose files
    # Assume there is only one distribution in the directory
    distribution_name = ""
    seeds::Vector{Int64} = []
    for f in files
        # Distribution name must end with [a-zA-Z]
        try
            s = split(split(f, "s=")[2], "_")[1]
            push!(seeds, parse(Int, s))
        catch
            continue
        end
    end
    if length(seeds)>0
        @logmsg settingLog "Found loose files, cleaning up... "
        clean_after_run(settings, seeds)
    end
end

function search_for_settings(path, dist)

    files = readdir(path*dist)
    settings = []
    for f in files

        # We skipp empty folders
        if isempty(readdir(path*dist*"/"*f))
            continue
        end
        params = split(f, " ")
        setting = Dict()

        for (key, value) in split.(params, "=")
            v_float = tryparse(Float64, value)
            v_int = tryparse(Int64, value)
            if v_int !== nothing
                setting[key] = v_int
            elseif v_float !== nothing
                setting[key] = v_float
            else
                setting[key] = value
            end
        end
        add_name_and_path_to_setting!(path, dist, setting)
        push!(settings, setting)
    end
    return settings
end

global_settings = nothing
function get_file_path(L, α, t, NR, dist="Uniform", data_path="data/"; average=true)
    setting = make_settings(L, t, NR, α, data_path, dist)
    return setting["path"]*setting["name"]*(average ? "" : "_bulk")*".jld2"
end



function load_file(L, α, t, NR, dist="Uniform"; data_path="data/", seed=-1, average=true)
    # We include this check so that we don't have to search for settings
    # every time we want to load a file
    if global_settings === nothing || data_path != "data/"
        global global_settings = search_for_settings(data_path, dist)
    end

    if NR=="LLS" || NR == "ELS"
        α=0.0
    end
    settings = filter(s -> s["L"] == L
                        && s["a"] == α 
                        && s["t"] == t
                        && s["nr"] == NR, global_settings)
    @assert length(settings) < 2 "There are multiple possibilities"
    if length(settings) == 0
        println("There is no file maching these settings α=$α nr=$NR, L=$L, t=$t")
        return nothing
    end
    setting = settings[1]
    return load(get_file_name(setting, seed, average))
end

function load_file(settings; seed=-1, average=true)
    s = settings
    path = split(s["path"], "/")[1]*"/"
    return load_file(s["L"], s["a"], s["t"], s["nr"], s["dist"],
                     data_path=path, seed=seed, average=average)
end

function remove_key(key, settings)
    #Does this even do anything?
    seeds = get_seeds_in_compact_file(get_name(settings))
    expand_file(settings)    
    filter!(e->e≠key,data_keys)
    condense_files(settings, seeds)
end

function add_key(key, value)
    # doesn't work
    settings = search_for_settings("data/", "Uniform")
    for average=[true, false], s = settings
        name = get_file_name(s, -1, average)
        f = load(name)
        if average
            if !haskey(f, "average_$key")
                f["average_$key"] = value
                JLD2.save(name, f)
            end
        else
            for s in f["seeds_used"]
                if !haskey(f, "$key/$s")
                    f["$key/$s"] = value
                    JLD2.save(name, f)
                end
            end
        end
    end
end


function rename(path)

    @logmsg settingLog "Searching for loose files... "
    files = readdir(path)
    # Find distribution with lose files
    # Assume there is only one distribution in the directory
    for f in files
        f_path = path*f*"/"        
        for ff in readdir("$f_path")
            ff_path = f_path*ff
            if occursin("a=2.0", ff) && occursin("nr=LLS", ff)
                mv(ff_path, f_path*replace(ff, "a=2.0" => "a=0.0"))
            end
        end
        
        if occursin("a=2.0", f) && occursin("nr=LLS", f)
            mv(f_path, replace(f_path, "a=2.0" => "a=0.0"))
        end
        
    end
end
#rename("data/Uniform/")
#add_key("simulation_time", 0)

function get_data_overview(path="data/", dists=["Uniform"])
    
    for dist in dists
        println("Data for: $dist")
        settings = search_for_settings(path, dist)
        s = []
        nr = []
        L = []
        t = []
        α = []
        for setting in settings
            seeds = get_seeds_in_file(setting)
            if isempty(seeds)
                seed_text = "none"
            else
                except = setdiff(minimum(seeds):maximum(seeds), seeds)
                seed_text = "$(minimum(seeds)) - $(maximum(seeds))" * (isempty(except) ? "" : " except $(join(except, " , "))")
            end
            push!(s, seed_text)
            push!(L, setting["L"])
            push!(nr, setting["nr"])
            push!(t, [setting["t"]])
            push!(α, setting["a"])
        end

        # Rearrange data
        data = [(nr[i], α[i], L[i], t[i], s[i]) for i in eachindex(s)]
        # Sort data
        sort!(data)
        # Gather t and s
        d = 0 #deleted elements
        for i in eachindex(data)
            i -= d
            if i>1 && data[i][5] == data[i-1][5] && data[i][3] == data[i-1][3]
                #if length(data[i-1][4]) == 2
                #    # That means that j should look something like 0.0 - 0.5
                #    data[i-1][4][2] = data[i][4][1]
                #else
                #end                
                push!(data[i-1][4], data[i][4][1])
                popat!(data, i)
                d += 1
            end
        end
        # Present data
        println("NR α  L  t  seeds")
        for i in eachindex(data)
            for j in eachindex(data[i])
                # j chooses the variable nr, α, L, t or s
                # If the value changes
                if i==1 || j>=4 || data[i][j] != data[i-1][j]
                    if j!=4
                        println("   "^(j-1)*"$(data[i][j])")
                    else
                        println("   "^(j-1)*"$(join(data[i][j], ", "))")
                    end
                end
            end
        end
    end
end

function get_bundle_from_file(file, L; nr="LLS", t=0.0, α=2.0, dist="Uniform", seed=1, progression=0, step=0,
                            without_storage=true, spanning=false, update_tension=true, return_simulation_time=false)
    if without_storage
        b = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist, without_storage=without_storage)
    else
        b,s = get_fb(L, seed, α=α, t=t, nr=nr, dist=dist, without_storage=without_storage)
        s.spanning_cluster_size_storage = file["spanning_cluster_size/$seed"]
        s.spanning_cluster_perimiter_storage = file["spanning_cluster_perimiter/$seed"]
        s.spanning_cluster_step = file["spanning_cluster_step/$seed"]
        s.most_stressed_fiber = file["most_stressed_fiber/$seed"]
        s.nr_clusters = file["nr_clusters/$seed"]
        s.largest_cluster = file["largest_cluster/$seed"]
        s.largest_perimiter = file["largest_perimiter/$seed"]
    end

    break_sequence = file["break_sequence/$seed"] 
    b.break_sequence[1:length(break_sequence)] = break_sequence
    simulation_time = file["simulation_time/$seed"]
    b.current_step= file["last_step/$seed"]
    if progression != 0
        @assert step==0
        break_sequence = break_sequence[1:round(Int, L*L*progression)]
        b.current_step = round(Int, L*L*progression)
    end
    if step > 0
        @assert progression==0
        b.current_step = step
        break_sequence = break_sequence[1:step]
    elseif step < 0
        @assert progression==0
        b.current_step = file["last_step/$seed"]+step
        break_sequence = break_sequence[1:(file["last_step/$seed"]+step)]
    end
    if spanning
        @assert step==0
        @assert progression==0
        b.current_step = file["spanning_cluster_step/$seed"]
        break_sequence = break_sequence[1:file["spanning_cluster_step/$seed"]]
    end

    break_fiber_list!(break_sequence, b)
    if update_tension
        update_tension!(b)
    else
        resetBundle!(b)
    end

    if without_storage
        return b
    else
        if return_simulation_time
            return b, s, simulation_time
        else
            return b, s
        end
    end
end

function get_bundles_from_settings(settings; seeds, progression=0, step=0,
        without_storage=true, update_tension=true, spanning=false, return_simulation_time=false)
    file = load_file(settings, average=false)
    L = settings["L"]
    nr = settings["nr"]
    t = settings["t"]
    α = settings["a"]
    dist = settings["dist"]
    N = L*L
    bundles = []
    for seed in seeds
        b = get_bundle_from_file(file, L, nr=nr, t=t, α=α, dist=dist, seed=seed, progression=progression,
            step=step, without_storage=without_storage, update_tension=update_tension, spanning=spanning, return_simulation_time=return_simulation_time)
        push!(bundles, b)
    end
    if length(seeds)==1
        return bundles[1]
    else
        return bundles
    end
end

function recalculate_average_file(path="data/", dists=["Uniform"]; max_seed=9999)
    # We want to expand all the files and then
    # condense all of them again
    # I found that it's not benificial to have many seeds, so we have the option to
    # delete seeds that are larger than max_seed
    return # Remove this line to actually use (This function is a bit dangerous)
    for dist in dists
        settings = search_for_settings(path, dist)
        @showprogress for s in settings
            search_for_loose_files(s)
            all_seeds = get_seeds_in_file(s)
            seeds = all_seeds[all_seeds .<= max_seed]
            unused_seeds = all_seeds[all_seeds .> max_seed]
            #println("$(s["L"]) ", length(all_seeds), ": ", length(seeds), ", ", length(unused_seeds))

            if isempty(seeds) || s["L"] > 16
                continue
            else
                expand_file(s)
                condense_files(s, seeds, keep_existing_seeds=false)
            end
        end
    end
    println("Success!")
end

function rename_files_and_folders(path="data/", dists=["gyration_data"])
    new_name(s) = replace(s, "r_slope" => "r")
    for dist in dists
        full_path = path*dist
        for file in readdir(full_path)
            new_file_name = new_name(file)
            if new_file_name != file
                mv(full_path*"/"*file, full_path*"/"*new_file_name)
            end
        end
    end
    println("Success!")
end

#rename_files_and_folders()

#recalculate_average_file()

#get_data_overview()
