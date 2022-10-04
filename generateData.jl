
include("dataGenerator.jl")

seeds = 1:11
L = [128]
t = (1:9) ./ 10
NHR = ["", "CNR", "SNR"]

itterate_settings(L, t, NHR, seeds; overwrite=false)