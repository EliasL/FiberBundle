using JLD2
using CodecLz4
using Logging
include("logingLevels.jl")

averaged_data_keys = [
                    "nr_clusters", 
                    "largest_cluster",
                    "largest_perimiter",
                    "most_stressed_fiber",
                    "spanning_cluster_size",
                    "spanning_cluster_perimiter",
                    "spanning_cluster_step",
                    ]
seed_specific_keys = [
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
        "α" => α,
    )
    settings["name"] = get_uniform_distribution_name(settings)
    settings["path"] = path*dist*"/"*settings["name"]*"/"

    if !isdir(settings["path"])
        mkpath(settings["path"])
    end
end

function get_uniform_distribution_name(settings)
    name = "Uniform"
    excluded = ["name", "path"]
    for (setting, value) in settings
        if setting in excluded
            continue
        else
            name *= " $setting=$value"
        end
    end
    return name
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

function get_seeds_in_compact_file(file_name)
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
    existing_seeds = get_seeds_in_compact_file(condensed_file_name)
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
                jldopen(seed_file_name, "r") do s_file
                    for key in averaged_data_keys
                        # If the key doesn't exist, make it and set it to zero
                        if !haskey(averages, key)
                            if length(s_file[key]) == 1
                                averages[key] = 0
                            else
                                averages[key] = zeros(length(s_file[key]))
                            end
                        end
                        averages[key] += s_file[key] ./ nr_seeds
                        condensed_file["$key/$seed"] = s_file[key]
                    end
                    if seed <= 10
                        for key in seed_specific_keys
                            condensed_file["$key/$seed"] = s_file[key]
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

    search_for_loose_files(path)


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

function search_for_loose_files(path)
    #TODO
    error("Not fixed yet")
    return
    @logmsg settingLog "Searching for loose files... "
    files = readdir(path)
    # Find distribution with lose files
    # Assume there is only one distribution in the directory
    distribution_name = ""
    Ls = Set([])
    seeds = Dict()
    for f in files
        # Distribution name must end with [a-zA-Z]
        m = match(r"(^.+[a-zA-Z])([0-9]+)+,([0-9]+)+_bulk.jld2$", f)
        if m !== nothing
            distribution_name = m.captures[1]
            L = parse(Int64, m.captures[2])
            seed = parse(Int64, m.captures[3])
            push!(Ls , L)
            if haskey(seeds, L)
                push!(seeds[L] , seed)
            else
                seeds[L] = [seed]
            end
        end
    end
    if length(Ls)>0
        @logmsg settingLog "Found loose files, cleaning up... "
    end
    for L in Ls
        clean_after_run(L, distribution_name, path, seeds[L])
    end
    @logmsg settingLog "Done!"
end

function search_for_t(path)

    files = readdir(path)
    t = Set([])
    for f in files
        # t must be on the form "t=number.number " ie. t=0.0 or t=12.34
        m = match(r"t=(([0-9]+)\.([0-9]+)) ", f)
        if m !== nothing
            push!(t, parse(Float64, m.captures[1]))
        end
    end
    return sort(collect(t))
end

function remove_key(key, settings)
    seeds = get_seeds_in_compact_file(get_name(settings))
    expand_file(settings)    
    filter!(e->e≠key,data_keys)
    condense_files(settings, seeds)
end 

function rename()

    @logmsg settingLog "Searching for loose files... "
    files = readdir("data")
    # Find distribution with lose files
    # Assume there is only one distribution in the directory
    for f in files
        t_unifor_nr = split(split(f, " ")
        t= parse(Float64, t_unifor_nr[1], "=")[2])
        nr = t_unifor_nr[3]
        
        for ff in readdir("data/$f")
            sett
            L_bulk = split(ff, nr)[2]
            L = split(L_bulk, "_")
            bulk = length(L) == length(L_bulk)


            settings = make_settings("Uniform", L, t, nr, 2, "data/")

            mv(ff, get_file_name(settings, average=!bulk))
            rm(ff)
        end
        rm(f)
    end

    distribution_name = ""
    Ls = Set([])
    seeds = Dict()
    for f in files
        # Distribution name must end with [a-zA-Z]
        m = match(r"(^.+[a-zA-Z])([0-9]+)+,([0-9]+)+_bulk.jld2$", f)
        if m !== nothing
            distribution_name = m.captures[1]
            L = parse(Int64, m.captures[2])
            seed = parse(Int64, m.captures[3])
            push!(Ls , L)
            if haskey(seeds, L)
                push!(seeds[L] , seed)
            else
                seeds[L] = [seed]
            end
        end
    end
    if length(Ls)>0
        @logmsg settingLog "Found loose files, cleaning up... "
    end
    for L in Ls
        clean_after_run(L, distribution_name, path, seeds[L])
    end
    @logmsg settingLog "Done!"
end