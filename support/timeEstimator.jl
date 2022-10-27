using Dates
using BenchmarkTools
using ProgressMeter
using Distributed

include("../dataGenerator.jl")
include("dataManager.jl")


include("logingLevels.jl")
# Sett logging level
logger = SimpleLogger(stdout, 0)
global_logger(logger)

function calculate_times()
    L = [16, 32, 64, 128, 256]
    NR = ["UNR", "SNR", "CNR"] # We will only use CNR since it takes the longest

    jldopen("data/time_estimates.jld2", "w") do file
        @showprogress for l=L, nr=NR
            # t is in seconds. Benchmark gives result in nano seconds, hence the 10^9
            file["$nr$l"] = mean(@benchmark run_setting($l, $nr)).time * 10^9 
        end
    end
end

function run_setting(L, NR)
    generate_data("", L, 1, "Uniform", 0.0, false, NR; save_data=false)
end

function time_estimate(dimensions, α, regimes, NRs, seeds, path)

    seconds_estimated = 0

    time_estimates = jldopen("data/time_estimates.jlds2", "r")

    for L=dimensions, t=regimes, nr=NRs, a=α
        distribution_name = get_uniform_distribution_name(t, nr)
        missing_seeds = get_missing_seeds(L, distribution_name, path, seeds)
        seconds_estimated += length(missing_seeds) * time_estimates["$nr$L"]
    end


    threads = Threads.nthreads()
    time_estimate = seconds_estimated/threads
    formated_time = Dates.canonicalize(Dates.CompoundPeriod(Dates.Second(floor(Int64, time_estimate))))
    @info "This will probably take: $formated_time"
end

calculate_times()
