using Distributed
using Suppressor: @suppress_err

# this just removes workers if there are leftovers from a crash
@suppress_err rmprocs(workers())

threads = Threads.nthreads()
print("Starting $threads workers... ")
addprocs(threads)
println("Done!")

include("dataGenerator.jl")

seeds = 0:500
L = [128, 64, 256]
t = vcat((0:8) ./ 20, (5:9) ./ 10)
NR = ["UNR", "CNR", "SNR"]

itterate_settings(L, t, NR, seeds; overwrite=true, estimate_time=true)

println("Removing workers")
rmprocs(workers())