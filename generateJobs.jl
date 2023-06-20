include("support/timeEstimator.jl")

function format(total_seconds)
    (d,r) = divrem(total_seconds,60*60*24)
    (h,r) = divrem(r,60*60)
    (m,r) = divrem(r, 60)
    s = ceil(Int64, r)
    (d,h,m,s) = round.(Int64, [d,h,m,s])
    extra_hours = floor(Int64, d*4 + h/6)
    extra_minutes = 15 + h*5
    extra_days=3

    return "$(d+extra_days)-$(h+extra_hours):$(m+extra_minutes):$s"
end

function make_job(s, L; t = (0:9) ./ 10, NR = ["CLS", "LLS"], α = [2.0], dist="Uniform", force_short=false)
    threads = 64
    seconds = time_estimate(L, α, t, NR, dist=dist, collect(seeds[1]:seeds[2]), threads=threads)
    formated_time = format(seconds)
    #partition = seconds < 3600 || force_short ? "short" : "porelab"
    if force_short
        formated_time = "0-1:0:0"
        partition = "short"
    else
        partition = "porelab"
end

    #SBATCH --exclusive=user
    file_text = """
    #!/bin/bash

    #SBATCH -J $(maximum(L))-$(seeds[2])
    #SBATCH -p $partition
    #SBATCH -N 1
    #SBATCH -n $threads
    #SBATCH --nice=50000
    #SBATCH --time=$formated_time

    ml eb
    ml Julia/1.7.2-linux-x86_64
    julia --threads $threads generateData.jl L $(join(L, " ")) t $(join(t, " ")) a $(join(α, " ")) NR $(join(NR, " ")) s $(join(seeds, " ")) dist $dist 

    wait
    """

    write("job.sh", file_text)
end

function start_job()
    run(`sbatch job.sh`)
end

t = vcat([0.5, 1],1.1:0.1:1.4, 1.5:0.5:5)
NR = ["LLS", "CLS"]
dist = "Weibull"

L=[64]
seeds = [0,2000]

#make_job(seeds, L, t=t, α=[2.0], NR=NR, dist=dist, force_short=false)
#start_job()

for single_t in t
    make_job(seeds, L, t=[single_t], α=[2.0], NR=NR, dist=dist)
    start_job()
end

#= L=[512]
seeds = [0,200]
make_job(seeds, L, t=[0.0], α=[2.0], NR=NR, dist=dist)
start_job() =#
# cancle jobs
#= for i in 1248451:1248483
    run(`scancel $i`)
end =#
