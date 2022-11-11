using Distributed
using Logging

include("support/logingLevels.jl")
# Sett logging level
logger = SimpleLogger(stdout, settingLog)
#logger = SimpleLogger(stdout, -10000)
global_logger(logger)





seeds = 0:500-1 # Zero indexing, -1 to get 1000 samples instead of 1001.
L = [128]
t = (0:0) ./ 10#vcat((0:8) ./ 20, (5:7) ./ 10, (16:19) ./20)
α = [2.0]#[1, 1.5, 2, 2.5, 3, 5, 9, 15, 30]
#t = vcat((0:8) ./ 20, (5:9) ./ 10)
NR = ["SNR", "UNR"]
use_threads = true
overwrite = false
#time_estimate(L, t, NR, seeds, rough_estimate=true)

if use_threads
    # this just removes workers if there are leftovers from a crash
    rmprocs(workers())

    threads = Threads.nthreads()
    @logmsg nodeLog "Starting $threads workers... "
    addprocs(threads; exeflags="--project=$(Base.active_project())")
    @everywhere include("dataGenerator.jl")

    @logmsg nodeLog "Running settings: \n $L, \n $t, \n $α"
    @time itterate_settings(L, α, t, NR, seeds; overwrite=overwrite)
    @logmsg nodeLog "Removing workers"
    rmprocs(workers())
else
    include("dataGenerator.jl")
    @logmsg nodeLog "Start run"
    itterate_settings(L, α, t, NR, seeds; overwrite=overwrite, use_threads=use_threads)

end