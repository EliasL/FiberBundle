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


function time_estimate(dimensions, α, regimes, NRs, seeds; path="data/", dist="Uniform", threads=0)

    seconds_estimated = 0

    for L=dimensions, t=regimes, nr=NRs, a=α
        settings_with_correct_a = make_settings(dist, L, t, nr, a, path)
        if a != 0.0
            # Any alpha should take roughly the same time
            a = 2.0
        end
        settings = make_settings(dist, L, t, nr, a, path)
        try
            missing_seeds = get_missing_seeds(settings_with_correct_a, seeds)
            f = load_file(settings)
            seconds_estimated += length(missing_seeds) * f["average_simulation_time"]
        catch
            @warn "Unable to estimate time, missing data for:\n$(display_setting(settings))"
            assume_1_day = 3600*24
            assume_7_days = 3600*24*7
            return assume_7_days
        end
    end

    if threads==0
        threads = Threads.nthreads()
    end
    if length(seeds) >= threads
        time_estimate = seconds_estimated/threads
    else
        time_estimate = seconds_estimated/length(seeds)
    end
    formated_time = Dates.canonicalize(Dates.CompoundPeriod(Dates.Second(floor(Int64, time_estimate))))
    if isempty(formated_time.periods)
        @info "This will probably not take very long."
    else
        @info "This will probably take: $formated_time"
    end
    return time_estimate
end

function test()
    seeds = 0:4000-1 # Zero indexing, -1 to get 1000 samples instead of 1001.
    L = [8, 16, 32, 64, 128]
    t = (0:9) ./ 10
    α = [2.0]#[1, 1.5, 2, 2.5, 3, 5, 9, 15, 30]
    #t = vcat((0:8) ./ 20, (5:9) ./ 10)
    NR = ["CLS", "LLS"]
    use_threads = true
    overwrite = false
    time_estimate(L, α, t, NR, seeds)
end
#test()