using JLD2

data_keys = ["nr_clusters", 
             "largest_cluster",
             "largest_perimiter",
             "most_stressed_fiber",
             "spanning_cluster_size",
             "spanning_cluster_perimiter",
             "spanning_cluster_step",
             # These are not averaged
             "sample_states",
             "tension",
             "spanning_cluster",]

# Put the keys you don't want averaged in here, and make sure they are at the back of the data_keys list
not_averaged_keys = ["sample_states","tension","spanning_cluster",] #TODO separate into to lists
nr_not_averaged_keys = length(not_averaged_keys)
nr_averaged_keys = length(data_keys) - nr_not_averaged_keys
averaged_keys = data_keys[1:nr_averaged_keys]
# These data_keys will not be averaged
seed_specific_keys = data_keys[nr_averaged_keys+1:end]

function get_name(L, distribution, path, seed::Int=-1, average=false)
    @assert seed>=-1 "We don't want to use negative seeds since -1 is special here"
    if average
        return path*distribution*"$L.jld2"
    elseif seed == -1
        return path*distribution*"$(L)_bulk.jld2"
    else
        return path*distribution*"$L,$(seed)_bulk.jld2"
    end
end

function make_get_name(L, distribution, path)
    return f(seed=-1; average=false) = get_name(L, distribution, path, seed, average)
end

function expand_file(L, distribution, path, overwritten_seeds::AbstractArray=Vector{Int64}([]))
    get_name_fun = make_get_name(L, distribution, path)
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

function condense_files(L, distribution, path, requested_seeds::AbstractArray; remove_files=true)
    get_name_fun = make_get_name(L, distribution, path)
    condensed_file_name = get_name_fun()
    averaged_file_name = get_name_fun(average=true)
    existing_seeds = get_seeds_in_compact_file(condensed_file_name)
    seeds = union(existing_seeds, requested_seeds)
    nr_seeds = length(seeds)

    averages = Dict()
    jldopen(condensed_file_name, "w") do condensed_file
    jldopen(averaged_file_name, "w") do averaged_file
        for seed in seeds
            seed_file_name = get_name_fun(seed)
            jldopen(seed_file_name, "r") do s_file
                for key in data_keys
                    if key in averaged_keys
                        # If the key doesn't exist, make it and set it to zero
                        if !haskey(averages, key)
                            if length(s_file[key]) == 1
                                averages[key] = 0
                            else
                                averages[key] = zeros(length(s_file[key]))
                            end
                        end
                        averages[key] += s_file[key] ./ nr_seeds
                    end
                    if key in not_averaged_keys && seed > 10
                        continue
                    end
                    condensed_file["$key/$seed"] = s_file[key]
                end
            end
            if remove_files
                rm(seed_file_name)
            end
        end
        for key in averaged_keys
            averaged_file["average_$key"] = averages[key]
        end
        averaged_file["seeds_used"] = seeds
        averaged_file["nr_seeds_used"] = length(seeds)
        condensed_file["seeds_used"] = seeds
        condensed_file["nr_seeds_used"] = length(seeds)
    end # Averaged file
    end # Condensed file
end

function get_missing_seeds(file_name, requested_seeds)

    if isfile(file_name)
        jldopen(file_name, "r") do existing_data
            existing_seeds = existing_data["seeds_used"]
            return Vector{Int64}(setdiff(requested_seeds, existing_seeds))
        end
    else
        return requested_seeds
    end
end

function prepare_run(L, distribution, path, requested_seeds::AbstractArray, overwrite=false)

    search_for_loose_files(path)


    get_name_fun = make_get_name(L, distribution, path)
    condensed_file_name = get_name_fun()
    # Get missing seeds
    if overwrite
        missing_seeds = requested_seeds
    else
        missing_seeds = get_missing_seeds(condensed_file_name, requested_seeds)
    end

    # Print some info
    nr_found_files = length(requested_seeds) - length(missing_seeds)
    if nr_found_files>0 && nr_found_files!=length(requested_seeds)
        println("$nr_found_files files were found and skipped.")
    end

    # Here we expand the existing datafile 
    if isfile(condensed_file_name)
        if length(missing_seeds) == 0
            println("All seeds already exist")
        else
            expand_file(L, distribution, path, missing_seeds)
        end
    end

    return missing_seeds
end

function clean_after_run(L, distribution, path, requested_seeds::AbstractArray)
    condense_files(L, distribution, path, requested_seeds, remove_files=true)
end

function search_for_loose_files(path)
    print("Searching for loose files... ")
    files = readdir(path)
    # Find distribution with lose files
    # Assume there is only one distribution in the directory
    distribution_name = ""
    Ls = Set([])
    seeds = Dict()
    for f in files
        # Distribution name must end with [a-zA-Z]
        m = match(r"(^.+[a-zA-Z])([0-9]+)+,([0-9]+)+.jld2$", f)
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
        print("\nFound loose files, cleaning up... ")
    end
    for L in Ls
        clean_after_run(L, distribution_name, path, seeds[L])
    end
    print("Done!\r")
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

function remove_key(key, L, distribution, path)
    seeds = get_seeds_in_compact_file(get_name(L, distribution, path))
    expand_file(L, distribution, path)    
    filter!(e->e≠key,data_keys)
    condense_files(L, distribution, path, seeds)
end 