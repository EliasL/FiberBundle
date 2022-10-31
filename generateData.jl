using Distributed
using Logging

include("support/logingLevels.jl")
# Sett logging level
logger = SimpleLogger(stdout, settingLog)
#logger = SimpleLogger(stdout, -10000)
global_logger(logger)





seeds = 0:48-1 # Zero indexing, -1 to get 1000 samples instead of 1001.
L = [32]
t = vcat((0:8) ./ 20, (5:7) ./ 10, (16:19) ./20)
α = 2.0#[1, 1.5, 2, 2.5, 3, 5, 9, 15, 30]
#t = vcat((0:8) ./ 20, (5:9) ./ 10)
NR = ["UNR", "CNR", "SNR"]
use_threads = true
#time_estimate(L, t, NR, seeds, rough_estimate=true)

if use_threads
    # this just removes workers if there are leftovers from a crash
    rmprocs(workers())

    threads = Threads.nthreads()
    @logmsg nodeLog "Starting $threads workers... "
    addprocs(threads; exeflags="--project=$(Base.active_project())")
    @logmsg nodeLog "Done!"

    @everywhere include("dataGenerator.jl")

    @logmsg nodeLog "Start run"
    @time itterate_settings(L, α, t, NR, seeds; overwrite=true)
    @logmsg nodeLog "Done"
    @logmsg nodeLog "Removing workers"
    rmprocs(workers())
else
    include("dataGenerator.jl")
    @logmsg nodeLog "Start run"
    itterate_settings(L, α, t, NR, seeds; overwrite=true, use_threads=use_threads)

end