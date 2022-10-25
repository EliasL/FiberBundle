using Distributed
using Suppressor: @suppress_err
using Logging
using MPI
MPI.Init()
comm = MPI.COMM_WORLD

include("support/logingLevels.jl")
# Sett logging level
logger = SimpleLogger(stdout, nodeLog)
global_logger(logger)





seeds = 0:200-1 # Zero indexing, -1 to get 1000 samples instead of 1001.
L = [64]
t = [0.1, 0.2]#vcat((0:8) ./ 20, (5:7) ./ 10, (16:19) ./20, [0.925, 0.975])
#t = vcat((0:8) ./ 20, (5:9) ./ 10)
NR = ["UNR"]#, "CNR", "SNR"]

#time_estimate(L, t, NR, seeds, rough_estimate=true)







print("Hello world, I am rank $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))\n")
if MPI.Comm_rank(comm) == 0
    @logmsg rootLog "I am groot"
    @logmsg rootLog "L=64 40 seeds 8 julia thread, 4 job threads"
end


# this just removes workers if there are leftovers from a crash
@suppress_err rmprocs(workers())

threads = Threads.nthreads()
@logmsg nodeLog "Starting $threads workers... "
addprocs(threads; exeflags="--project=$(Base.active_project())")
@logmsg nodeLog "Done!"

include("dataGenerator.jl")

@logmsg nodeLog "Start run"
@time itterate_settings(L, t, NR, seeds; overwrite=true)
@logmsg nodeLog "Done"
@logmsg nodeLog "Removing workers"
rmprocs(workers())

MPI.Barrier(comm)