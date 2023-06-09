using LinearAlgebra
using Interpolations

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

function σ_to_energy(σ_max)
    @warn "Not sure if works! Check dddΔ"
    N = length(σ_max)
    d = (1:N) ./ N #Damage
    κ = 1 # Spring constant
    # NB! This function assumes that σ_max is equal to stretch. For that to be true,
    # κ is required to be 1!
    Δ = σ_max # Stretch/elongation
    p = 1
    #dddΔ =  #derivative of damage with respect to elongation
    #E = @. N*κ/2 *(1 -d)*Δ^2 
    E = @. N*κ/2 * (2*Δ*(1-d) - Δ^2*p)
    dEdΔ = N*κ .*Δ .*(1 .-d)
    return E, dEdΔ
end

function find_localization(cluster_size, critical_gradiant=1)
    # gives the first x at which the cluster is growing at a certain rate
    interp = interpolate((eachindex(cluster_size),), cluster_size*length(cluster_size), Gridded(Linear()))
    grad = only.(Interpolations.gradient.(Ref(interp), eachindex(cluster_size)))[10:end]
    #p = plot(grad)
    #display(p)
    #critical_gradiant = (maximum(grad) - minimum(grad))/2
    return findfirst(x->x>=critical_gradiant, grad)
end

function find_localization_nr_clusters(nr_clusters)
    @warn "This does not work well!" #It does find the max of the number of clusters, but
    #= # that does not corespond to localization (It seems)
    # gives the first x at which the number of clusters stop growing
    interp = interpolate((eachindex(nr_clusters),), nr_clusters, Gridded(Linear()))
    grad = only.(Interpolations.gradient.(Ref(interp), eachindex(nr_clusters)))[10:end]
    #p = plot(grad)
    #display(p)
    critical_gradiant = 0
    return findfirst(x-> isapprox(x,critical_gradiant;atol=10^-7), grad)  =#
    return argmax(nr_clusters)
end

function find_avalanches(σ)
    # σ is the most stressed fiber, so max(σ), not all N of them, but one for each k
    #TODO test function
    current_avalanche_index = 0
    in_avalanche =  false    
    a = zeros(Int64, length(σ))
    for i in eachindex(σ)[1:end-1]
        if σ[i] > σ[i+1] && in_avalanche
            current_avalanche_index=i
            in_avalanche = true
        else
            a[i] = 0
        end
        if σ[i] < σ[i+1] && in_avalanvhe
            a[current_avalanche_index] = i - current_avalanche_index
            in_avalanche = false
        end
    end
    return a
end


function test()
    nr = "ELS"
    path = "data/"
    t = 0.0
    L=128
    α = 2.0
    seed = 1
    settings = make_settings(L, t, nr, α)
    b = get_bundles_from_settings(settings, seeds=seed, spanning=true)
    p = plot_fb(b, show=false)
    minor_axes, major_axes, minor_values, major_values = find_major_and_minor_axes(b)
    plot_fb_axes(b, minor_axes, major_axes, minor_values, major_values)
    R = find_radius_of_gyration(b)
    plot_gyration_radi(b, R, nr=1)
    display(p)
end
#test()

function box_counting(b::FB)
    # If L=32, returns conunts for boxes of size 2, 4, 8 and 16, but not 1 and 32
    # Find spanning cluster
    @assert b.spanning_cluster_id != -1
    sc = b.spanning_cluster_id

    # box counts
    nr_sizes = round(Int64, log2(b.L))-1

    s_boxes = []
    h_boxes = []
    # Here we create all the boxes. Later we count how many are true
    for boxes in [s_boxes, h_boxes]
        for i in 1:nr_sizes
            push!(boxes, zeros(Bool, (2^(nr_sizes - i+1), 2^(nr_sizes - i+1))))
        end
    end

    for (i, state) in enumerate(b.status)

        if state == sc
            x, y = fiber_index_to_xy(i, b.L)
            #println("$x, $y")
            for i in 1:nr_sizes
                pow = 2^i
                # We want to find the coordinate for each of the larger boxes
                new_x = ceil(Int64, x/pow)
                new_y = ceil(Int64, y/pow)
                s_boxes[i][new_x, new_y] = true
            end
        elseif sc in b.status[b.neighbours[i, :]] # state != sc
            x, y = fiber_index_to_xy(i, b.L)
            #println("$x, $y")
            for i in 1:nr_sizes
                pow = 2^i
                # We want to find the coordinate for each of the larger boxes
                new_x = ceil(Int64, x/pow)
                new_y = ceil(Int64, y/pow)
                h_boxes[i][new_x, new_y] = true
            end
        end

    end


    #plot_fb(b, use_shift=false)
    #display.(heatmap.(rotl90.(s_boxes)))
    #display.(heatmap.(rotl90.(h_boxes)))
    s_count = [count(a) for a in s_boxes]
    h_count = [count(a) for a in h_boxes]
    pushfirst!(s_count, b.cluster_size[sc])
    pushfirst!(h_count, b.cluster_outline_length[sc])
    return s_count, h_count
end