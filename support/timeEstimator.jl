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


function time_estimate(dimensions, α, regimes, NRs, seeds, path="data/", dist="Uniform", threads=0)

    seconds_estimated = 0

    for L=dimensions, t=regimes, nr=NRs, a=α
        settings = make_settings(dist, L, t, nr, a, path)
        try
            missing_seeds = get_missing_seeds(settings, seeds)
            f = load_file(settings)
            seconds_estimated += length(missing_seeds) * f["average_simulation_time"]
        catch
            @warn "Unable to estimate time, missing data for $settings"
            return
        end
    end

    if threads==0
        threads = Threads.nthreads()
    end
    time_estimate = seconds_estimated/threads

    formated_time = Dates.canonicalize(Dates.CompoundPeriod(Dates.Second(floor(Int64, time_estimate))))
    @info "This will probably take: $formated_time"
end
