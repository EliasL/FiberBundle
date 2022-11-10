using LinearAlgebra

include("dataManager.jl")
include("../burningMan.jl")
include("../plottingScripts/showBundle.jl")


function inertia_tensors_of_clusters(b::FB)
    I = zeros(Float64, b.c, 3)
    for (i, state) in enumerate(b.status)
        if state > BROKEN # Then the fiber is part of a cluser
            x, y = fiber_index_to_xy(i, b.L)
            cmx = b.cluster_cm_x[state]
            cmy = b.cluster_cm_y[state]
            #Ixx
            I[state, 1] += (cmy - y)^2
            #Iyy
            I[state, 2] += (cmx - x)^2
            #Ixy
            I[state, 3] -= (cmy - y) * (cmy - x)
        end
    end
    return I
end

function diagonalize_inertia_tensors(I_::AbstractMatrix{Float64})
    nr_clusters = size(I_)[1]
    minor_axes = zeros(Float64, nr_clusters, 2)
    major_axes = zeros(Float64, nr_clusters, 2)

    for (i, I) in enumerate(eachrow(I_))
        I = SymTridiagonal(I[1:2], I[3:3])
        d = eigen(I)
        if d.values[1] > d.values[2]
            major_axes[i] = d.vectors[1]
            minor_axes[i] = d.vectors[2]
        else
            major_axes[i] = d.vectors[2]
            minor_axes[i] = d.vectors[1]
        end
    end
    return minor_axes, major_axes
end 

function find_major_and_minor_axes(b::FB)
    I = inertia_tensors_of_clusters(b)
    return diagonalize_inertia_tensors(I)
end

function test()
    nr = "SNR"
    path = "data/"
    t = 0.1
    L=256
    α = 2.0
    seed = 1
    setting = make_settings("Uniform", L, t, nr, α, path)
    file = load_file(setting, average=false)
    b = get_fb(L, nr=nr, without_storage=true)
    b.status = file["spanning_cluster_state/$seed"]
    shift_spanning_cluster!(b)
    resetBundle!(b)
    update_σ!(b)
    minor_axes, major_axes = find_major_and_minor_axies(b)
    plot_fb(b, show=false)
end

test()