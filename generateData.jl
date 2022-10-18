using Distributed
using Suppressor: @suppress_err

# this just removes workers if there are leftovers from a crash
@suppress_err rmprocs(workers())

threads = Threads.nthreads()
print("Starting $threads workers... ")
addprocs(threads; exeflags="--project=$(Base.active_project())")
println("Done!")

include("dataGenerator.jl")

seeds = 0:10-1 # Zero indexing, -1 to get 1000 samples instead of 1001.
L = [32]
t = vcat((0:8) ./ 20, (5:7) ./ 10, (16:19) ./20, [0.925, 0.975])
NR = ["UNR", "CNR", "SNR"]

#time_estimate(L, t, NR, seeds, rough_estimate=true)

itterate_settings(L, t, NR, seeds; overwrite=false, estimate_time=false)

println("Removing workers")
rmprocs(workers())