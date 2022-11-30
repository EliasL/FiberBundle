include("support/timeEstimator.jl")

function format(total_seconds)
    (d,r) = divrem(total_seconds,60*60*24)
    (h,r) = divrem(r,60*60)
    (m,r) = divrem(r, 60)
    s = ceil(Int64, r)
    (d,h,m,s) = round.(Int64, [d,h,m,s])
    extra_hours = floor(Int64, d*2 + h/10)
    extra_minutes = 15 

    return "$d-$h:$m:$s"
end

function make_job(s, L; t=t = (0:9) ./ 10, NR = ["SNR", "UNR"], α = [2.0], force_short=false)
    threads = 40
    seconds = time_estimate(L, α, t, NR, collect(seeds[1]:seeds[2]), threads=threads)
    formated_time = format(seconds)
    partition = seconds < 3600 || force_short ? "short" : "porelab"
    if force_short
        formated_time = "0-1:0:0"
    end

    file_text = """
    #!/bin/bash

    #SBATCH -J $(maximum(L))-$(seeds[2])
    #SBATCH -p $partition
    #SBATCH -N 1
    #SBATCH -n 64
    #SBATCH --exclusive=user
    #SBATCH --time=$formated_time

    ml eb
    ml Julia/1.7.2-linux-x86_64
    julia --threads $threads generateData.jl L $(join(L, " ")) t $(join(t, " ")) NR $(join(NR, " ")) s $(join(seeds, " ")) 

    wait
    """

    write("job.sh", file_text)
end

function start_job()
    run(`sbatch job.sh`)
end



seeds = [0, 4000] # From seed to seed
L = [8,16,32,64,128]

make_job(seeds, L, force_short=false)
start_job()