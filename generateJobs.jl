include("support/timeEstimator.jl")

function format(date::Dates.CompoundPeriod)
    
    @assert length(date.periods) <= 4 "This will take over a month..."
    int_times = map(x -> x.value, date.periods)
    t = vcat(zeros(Int64, 4-length(int_times)), int_times)
    extra_hours = 1
    return "$(t[1])-$(t[2]+extra_hours):$(t[3]):$(t[4])"
end

function make_job(s, L, t=t = (0:9) ./ 10, NR = ["SNR", "UNR"], α = [2.0])
    threads = 40
    time = time_estimate(L, α, t, NR, collect(seeds[1]:seeds[2]), threads=threads)
    formated_time = format(time)
    file_text = """
    #!/bin/bash

    #SBATCH -J $(maximum(L)), $(seeds[2])
    #SBATCH -p porelab
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



seeds = [0, 100] # From seed to seed
L = [8, 16, 32, 64, 128]

make_job(seeds, L)
start_job()