using Distributed
using Logging

include("support/logingLevels.jl")
include("support/timeEstimator.jl")
# Sett logging level
logger = SimpleLogger(stdout, settingLog)
#logger = SimpleLogger(stdout, -10000)
global_logger(logger)

function get_ARGS()
    args = Dict()
    args["L"]=[]
    args["t"]=[]
    args["NR"]=[]
    args["s"]=[]
    current_setting = nothing
    for value in ARGS
        if value in keys(args)
            current_setting = value
        else
            if current_setting in ["s", "L"]
                push!(args[current_setting], parse(Int64, value))
            elseif current_setting in ["t"]
                push!(args[current_setting], parse(Float64, value))
            else
                @assert current_setting !== nothing "Invalid argument. Not one of $(keys(args))"
                push!(args[current_setting], value)
            end
        end
    end
    return args
end

#seeds = 0:100-1 # Zero indexing, -1 to get 1000 samples instead of 1001.
#L = [256, 512]
#t = (0:9) ./ 10
#NR = ["SNR", "UNR"]
args = get_ARGS()
L=args["L"]
t=args["t"]
NR=args["NR"]
s=args["s"]
s = s[1]:(s[2]-1) # Zero indexing, -1 to get 1000 samples instead of 1001.
α = [2.0]
use_threads = true
overwrite = false
#time_estimate(L, α, t, NR, seeds)

if use_threads
    # this just removes workers if there are leftovers from a crash
    rmprocs(workers())

    threads = Threads.nthreads()
    @logmsg nodeLog "Starting $threads workers... "
    addprocs(threads; exeflags="--project=$(Base.active_project())")
    @everywhere include("dataGenerator.jl")

    @logmsg nodeLog "Running settings: \n $L, \n $t, \n $α"
    @time itterate_settings(L, α, t, NR, s; overwrite=overwrite)
    @logmsg nodeLog "Removing workers"
    rmprocs(workers())
else
    include("dataGenerator.jl")
    @logmsg nodeLog "Start run"
    itterate_settings(L, α, t, NR, seeds; overwrite=overwrite, use_threads=use_threads)
end