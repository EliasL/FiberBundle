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
    threads = 50
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
    #SBATCH -n 64
    #SBATCH --nice
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



seeds = [0, 1000] # From seed to seed
L = [256]
#t = [0.38]
#t = vcat((0:9) ./ 10)
#t = vcat((0:10) ./ 50)
#t = vcat((11:20) ./ 50)
t = vcat((0:9) ./ 50, (10:20) ./ 50, (5:9) ./ 10)
#t = vcat((0:5) ./ 50)
#t = vcat((6:10) ./ 50)
#t = vcat((11:15) ./ 50)
#t = vcat((16:20) ./ 50)
#NB Alpha in code should be one higher than in the paper! α=2 in code means α=1 in paper.
dist = "ConstantAverageUniform"
#seeds = [0, 100] # From seed to seed
#make_job(seeds, L, t=t, α=[2.0], NR=NR, dist=dist, force_short=false)
#start_job()


#L = [64]
#make_job(seeds, L, t=t, α=[2.0], NR=NR, force_short=false)
#start_job()

NR = ["LLS", "CLS"]
for ts in t
    make_job(seeds, L, t=[ts], α=[2.0], NR=NR, dist=dist, force_short=false)
    start_job()
end
#t = vcat((6:10) ./ 50)

#make_job(seeds, L, t=t, α=[2.0], NR=NR, dist=dist, force_short=false)
#start_job()


#=
t = vcat((11:15) ./ 50)
make_job(seeds, L, t=t, α=[2.0], NR=NR, force_short=false)
start_job()

t = vcat((16:20) ./ 50)
make_job(seeds, L, t=t, α=[2.0], NR=NR, force_short=false)
start_job()
 =#