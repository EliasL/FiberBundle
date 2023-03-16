include("dataManager.jl")


function get_data_kN(L, NR, ts, dist, key; average=true, divide=:N, return_kN=true, data_path="newData")
    data = []
    kN = []
    for l in L
        N = l*l
        if divide == :N
            divisor=N
        else
            divisor=divide
        end
        LData = zeros(N, length(ts), length(NR))
        LkN = zeros(size(LData))
        for nr=eachindex(NR), t=eachindex(ts)
            name = get_file_name(l, 2.0, ts[t], NR[nr], dist, average=average, data_path=data_path)
            jldopen(name, "r") do file
                fdata = file[key]
                LData[:, t, nr] .= fdata./divisor
                LkN[:, t, nr] .= (1:N)./N
            end
        end          
        push!(data, LData)
        push!(kN, LkN)
    end
    if return_kN
        return data, kN
    else
        return data
    end
end


function get_data(L, NR, ts, dist, key, xKey, xKeyf; average=true, ex=[2,2], α=2.0, return_x=false, data_path="newData")
    do_data_test = true
    data = zeros(length(ts), length(L), length(NR))
    x_values = zeros(length(ts), length(L), length(NR))
    for t=eachindex(ts), l=eachindex(L), nr=eachindex(NR)
        if do_data_test
            data_test(L[l], α, ts[t], NR[nr], dist)
        end
        name = get_file_name(L[l], α, ts[t], NR[nr], dist, average=average, data_path=data_path)
            
        jldopen(name, "r") do file            
            if average==false
                value, avg_x = load_data(file, key, L[l], NR[nr], ts[t], dist, xKey, xKeyf, data_path)
                x = avg_x
            else
                x = xKeyf(file["average_$xKey"])
                value = file[key][round(Int64, x)]
            end
            data[t, l, nr] = value/L[l]^ex[nr]
            x_values[t, l, nr] = x
        end
    end        
    if return_x
        return data, x_values
    else 
        return data
    end
end


function load_data(bulk_file, key, l, nr, t, dist, xKey, xKeyf,data_path)
    path = "$data_path/$xKey/"
    name = "$(dist)_$(l)_$(nr)_$(t)_$key"
    avg_x = 0
    if !isdir(path)
        mkpath(path)
    end
    if isfile("$(path)$(name).jld2")
        f = load("$(path)$(name).jld2")
        value = f[name]
        avg_x = f[name*"x"]
    else
        value = 0
        seeds_used = bulk_file["seeds_used"]
        nr_seeds = length(seeds_used)
        for seed in bulk_file["seeds_used"]
            x = xKeyf(bulk_file["$xKey/$seed"])
            value += bulk_file["$key/$seed"][x]
            avg_x += x
        end
        value /= nr_seeds
        avg_x /= nr_seeds
        # Save data
        jldopen("$(path)$(name).jld2", "w") do file
            file[name] = value
            file[name*"x"] = avg_x
        end
    end

    return value, avg_x
end


function data_test(L, α, t, NR, dist, data_path="newData/")
    name = get_file_name(L, α, t, NR, dist, average=false, data_path=data_path)
    jldopen(name, "r") do file
        nr_seeds = file["nr_seeds_used"]
        if file["last_step/$(nr_seeds-1)"] != L*L
            @warn "Not fully broken: L=$L, t0=$t, NR=$NR"
        end
    end
end
#d = get_data([128], ["CLS"], [0.0, 0.1], [[1,2]], "average_most_stressed_fiber")
#println(d)