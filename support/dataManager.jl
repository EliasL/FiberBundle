using JLD2
using CodecLz4
using Logging
include("logingLevels.jl")

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
                    "break_sequence",
                    "sample_states",
                    "tension",
                    "spanning_cluster_state",
                    "spanning_cluster_tension",
                    ]
data_keys = vcat(averaged_data_keys, seed_specific_keys)


function make_settings(dist, L, t, nr, α, path)
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

function expand_file(settings, overwritten_seeds::AbstractArray=Vector{Int64}([]))
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
            # Write out the file
            jldopen(get_name_fun(seed), "w") do s_file
                for key in data_keys
                    if haskey(file, "$key/$seed")
                        s_file[key] = file["$key/$seed"]
                    end
                end
            end # Close seed file
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

function condense_files(settings, requested_seeds::AbstractArray; remove_files=true)
    get_name_fun = make_get_name(settings)
    condensed_file_name = get_name_fun()
    averaged_file_name = get_name_fun(average=true)
    existing_seeds = get_seeds_in_file(settings, )
    if length(existing_seeds) > 0
        expand_file(settings)
    end
    seeds = union(existing_seeds, requested_seeds)
    nr_seeds = length(seeds)

    averages = Dict()
    #try
        jldopen(condensed_file_name*".temp", "w", compress = LZ4FrameCompressor()) do condensed_file
        jldopen(averaged_file_name*".temp", "w") do averaged_file
            averaged_file["seeds_used"] = seeds
            averaged_file["nr_seeds_used"] = length(seeds)
            condensed_file["seeds_used"] = seeds
            condensed_file["nr_seeds_used"] = length(seeds)
            for seed in seeds
                seed_file_name = get_name_fun(seed)
                jldopen(seed_file_name, "r+") do s_file
                    for key in averaged_data_keys
                        if !haskey(s_file, key)
                            value = 0
                            @warn "$key not found in $(seed_file_name)!"
                        else
                            value = s_file[key]
                        end
                        # If the key doesn't exist, make it and set it to zero
                        if !haskey(averages, key)
                            if length(value) == 1
                                averages[key] = 0
                            else
                                averages[key] = zeros(length(value))
                            end
                        end
                        #println(value)
                        #println(averages[key])
                        averages[key] += value ./ nr_seeds
                        condensed_file["$key/$seed"] = value
                    end
                    if seed <= 10
                        for key in seed_specific_keys
                            if !haskey(s_file, key)
                                value = 0
                                @warn "$key not found in $(seed_file_name)!"
                            else
                                value = s_file[key]
                            end
                            condensed_file["$key/$seed"] = value
                        end
                    end
                end
            end
            for key in averaged_data_keys
                averaged_file["average_$key"] = averages[key]
                #TODO add uncertainty
            end
        end # Averaged file
        end # Condensed file
        if remove_files
            for seed in seeds
                rm(get_name_fun(seed))
            end
        end
        mv(condensed_file_name*".temp", condensed_file_name, force=true)
        mv(averaged_file_name*".temp", averaged_file_name, force=true)
    #catch e
    #    rm(condensed_file_name*".temp")
    #    rm(averaged_file_name*".temp")
    #    println("File condensing failed!")
    #    if e != "FileNot found"
    #        throw(e)
    #    end
    #    println(e)
    #end
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
    seeds = []
    for f in files
        # Distribution name must end with [a-zA-Z]
        try
            s = split(split(f, "s=")[2], "_")[1]
            push!(seeds, s)
        catch
            continue
        end
    end
    if length(seeds)>0
        @logmsg settingLog "Found loose files, cleaning up... "
        clean_after_run(settings, seeds)
    end
    @logmsg settingLog "Done!"
end

function search_for_settings(path, dist)

    files = readdir(path*dist)
    settings = []
    for f in files
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
    setting = make_settings(dist, L, t, NR, α, data_path)
    return setting["path"]*setting["name"]*(average ? "" : "_bulk")*".jld2"
end

function load_file(L, α, t, NR, dist="Uniform"; data_path="data/", seed=-1, average=true)
    # We include this check so that we don't have to search for settings
    # every time we want to load a file
    if global_settings === nothing
        global global_settings = search_for_settings(data_path, dist)
    end

    if NR=="UNR"
        α=0.0
    end
    settings = filter(s -> s["L"] == L
                        && s["a"]==α 
                        && s["t"]==t
                        && s["nr"]==NR, global_settings)
    @assert length(settings) < 2 "There are multiple possibilities"
    @assert length(settings) != 0 "There is no file maching these settings α=$α nr=$NR, L=$L, t=$t"
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
        println(name)
        
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
            if occursin("a=2.0", ff) && occursin("nr=UNR", ff)
                mv(ff_path, f_path*replace(ff, "a=2.0" => "a=0.0"))
            end
        end
        
        if occursin("a=2.0", f) && occursin("nr=UNR", f)
            mv(f_path, replace(f_path, "a=2.0" => "a=0.0"))
        end
        
    end
end
#rename("data/Uniform/")
#add_key("simulation_time", 0)

function get_data_overview(path="data/", dists=["Uniform"])
    nr = []
    L = []
    t = []
    α = []

    for dist in dists
        println("Data for: $dist")
        settings = search_for_settings(path, dist)
        for setting in settings
            s = get_seeds_in_file(setting)
            L = setting["L"]
            nr = setting["nr"]
            t = setting["t"]
            α = setting["a"]
            
            if length(s) == 0
                println("Empty: $L, $nr, $t, $α")
                continue
            end
            except = setdiff(minimum(s):maximum(s), s)
            seed_text = "$(minimum(s)) - $(maximum(s))" * (isempty(except) ? "" : " except $(join(except, " , "))")

            println("$L, $nr, $t, $α, $seed_text")
        end
    end
end

#get_data_overview()
