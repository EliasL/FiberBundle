include("dataManager.jl")


function get_data_kN(L, NR, ts, key; average=true, divide=:N, return_kN=true)
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
            name = get_file_name(l, 2.0, ts[t], NR[nr], "Uniform", average=average)
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



function get_data_spanning(L, NR, ts, key; average=true, divide=:N)
    do_data_test = true
    data = zeros(length(ts), length(L), length(NR))
    for t=eachindex(ts), l=eachindex(L), nr=eachindex(NR)
        if do_data_test
            data_test(L[l], 2.0, ts[t], NR[nr], "Uniform")
        end
        name = get_file_name(L[l], 2.0, ts[t], NR[nr], "Uniform", average=average)
        jldopen(name, "r") do file

            
            N = L[l]^2
            x = file["average_spanning_cluster_step"]
            #x = N*0.4
            d = file[key][round(Int64, x)]
            if divide == :N
                divisor=N
            elseif divide == :L
                divisor = L[l]
            elseif divide in [:max, :min, :maxSubMin]
                divisor = 1
            else
                divisor=divide
            end
            data[t, l, nr] = d/divisor
            if t == length(ts)
                if divide in [:max, :min, :maxSubMin]
                    f = divide == :min ? minimum : maximum 
                    # When we are at the last t value, we find the max of
                    # all the t values and normalize
                    max_t = f(data[:, l, nr])
                    if divide == :maxSubMin
                        data[:, l, nr] .-= minimum(data[:, l, nr])
                    end
                    data[:, l, nr] ./= max_t
                end
            end
        end
    end        
    return data
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