include("dataManager.jl")


function get_data_kN(L, NR, ts, dist, key; average=true, divide=:N, return_kN=true)
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
            name = get_file_name(l, 2.0, ts[t], NR[nr], dist, average=average)
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


function get_data(L, NR, ts, dist, key, xKey, xKeyf; average=true, ex=[2,2], α=2.0)
    do_data_test = true
    data = zeros(length(ts), length(L), length(NR))
    for t=eachindex(ts), l=eachindex(L), nr=eachindex(NR)
        if do_data_test
            data_test(L[l], α, ts[t], NR[nr], dist)
        end
        name = get_file_name(L[l], α, ts[t], NR[nr], dist, average=average)
        jldopen(name, "r") do file            
            if average==false
                value = load_data(file, key, L[l], NR[nr], ts[t], xKey, xKeyf)
            else
                x = xKeyf(file["average_$xKey"])
                value = file[key][ceil(Int64, x)]
            end
            data[t, l, nr] = value/L[l]^ex[nr]
        end
    end        
    return data
end



function load_data(bulk_file, key, l, nr, t, xKey, xKeyf)
    path = "data/$xKey/"
    if !isdir(path)
        mkpath(path)
    end
    if isfile("$(path)$(l)_$(nr)_$(t)_$key.jld2")
        f = load("$(path)$(l)_$(nr)_$(t)_$key.jld2")
        value = f["$(l)_$(nr)_$(t)_$key"]
    else
        value = 0
        seeds_used = bulk_file["seeds_used"]
        nr_seeds = length(seeds_used)
        for seed in bulk_file["seeds_used"]
            x = xKeyf(bulk_file["$xKey/$seed"])
            value += bulk_file["$key/$seed"][x]
        end
        value /= nr_seeds
        # Save data
        jldopen("$(path)$(l)_$(nr)_$(t)_$key.jld2", "w") do file
            file["$(l)_$(nr)_$(t)_$key"] = value
        end
    end

    return value
end


function data_test(L, α, t, NR, dist)
    name = get_file_name(L, α, t, NR, dist, average=false)
    jldopen(name, "r") do file
        nr_seeds = file["nr_seeds_used"]
        if file["last_step/$(nr_seeds-1)"] != L*L
            @warn "Not fully broken: L=$L, t0=$t, NR=$NR"
        end
    end
end
#d = get_data([128], ["CLS"], [0.0, 0.1], [[1,2]], "average_most_stressed_fiber")
#println(d)