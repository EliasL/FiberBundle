using LinearAlgebra

include("dataManager.jl")
include("../burningMan.jl")
include("../plottingScripts/showBundle.jl")

d = distance

function inertia_tensors_of_clusters(b::FB)
    I = zeros(Float64, b.c, 3)
    for (i, state) in enumerate(b.status)
        if state > 0#BROKEN # Then the fiber is part of a cluser
            x, y = fiber_index_to_xy(i, b.L)
            cmx = b.cluster_cm_x[state]
            cmy = b.cluster_cm_y[state]
            #Ixx
            I[state, 1] += (d(cmy, y, b.L))^2
            #Iyy
            I[state, 2] += (d(cmx, x, b.L))^2
            #Ixy
            I[state, 3] -= d(cmy, y, b.L) * d(cmx, x, b.L)
        end
    end
    return I
end

function diagonalize_inertia_tensors(I_::AbstractMatrix{Float64})
    nr_clusters = size(I_)[1]
    minor_axes = zeros(Float64, nr_clusters, 2)
    major_axes = zeros(Float64, nr_clusters, 2)
    minor_values = zeros(Float64, nr_clusters)
    major_values = zeros(Float64, nr_clusters)

    for (i, I) in enumerate(eachrow(I_))
        I = SymTridiagonal(I[1:2], I[3:3])
        d = eigen(I)
        if d.values[1] > d.values[2]
            major_axes[i,:] = d.vectors[1,:]
            major_values[i] = d.values[1]
            minor_axes[i,:] = d.vectors[2,:]
            minor_values[i] = d.values[2]
        else
            minor_axes[i,:] = d.vectors[1,:]
            minor_values[i] = d.values[1] 
            major_axes[i,:] = d.vectors[2,:]
            major_values[i] = d.values[2]
        end
    end
    return minor_axes, major_axes, minor_values, major_values
end 

function find_major_and_minor_axes(b::FB)
    I = inertia_tensors_of_clusters(b)
    return diagonalize_inertia_tensors(I)
end

function find_radius_of_gyration(b::FB)
    R = zeros(Float64, b.c)
    for (i, state) in enumerate(b.status)
        if state > 0#BROKEN # Then the fiber is part of a cluser
            x, y = fiber_index_to_xy(i, b.L)
            cmx = b.cluster_cm_x[state]
            cmy = b.cluster_cm_y[state]
            R[state] += d(cmy, y, b.L)^2 + d(cmx, x, b.L)^2
        end
    end
    R .= sqrt.(R ./ view(b.cluster_size, 1:b.c)) 
    return R
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
    minor_axes, major_axes, minor_values, major_values = find_major_and_minor_axes(b)
    p = plot_fb(b, show=false)
    #plot_fb_axes(b, minor_axes, major_axes, minor_values, major_values)
    R = find_radius_of_gyration(b)
    plot_gyration_radi(b, R)
    display(p)
end

test()