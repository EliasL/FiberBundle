
include("dataGenerator.jl")

seeds = 0:30
L = [64]
t = (0:9) ./ 10
NR = ["UNR", "CNR", "SNR"]

itterate_settings(L, t, NR, seeds; overwrite=true)