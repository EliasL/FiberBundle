using StaticArrays

Base.@kwdef struct neighbourhoodData{ A, B, F<:AbstractFloat, I<:Integer}
    L::I
    N::F
    x::MVector{A,I}
    y::MVector{B,I}
end

#test = neighbourhoodData(0x1, 0.0, [0x1,0x0], [0x1,0x0, 0x3])
a = SVector{3,typeof(1)}([1,2,3])
b = SVector{4,typeof(1)}([1,2,3, 4])
#test = neighbourhoodData{3,4}(1, 0.0, a, b)

Base.@kwdef struct Foo3{l, n, A, F<:AbstractFloat, I<:Integer}
    L::I = l
    N::I = n
    t::I = 4
    alpha::SVector{n,F} = SVector{n}(zeros(F, n))
    beta::SMatrix{n,A,I} = SMatrix{n,A}(zeros(I,n, A))
    beta2::SMatrix{n,2,I} = SMatrix{n,2}(zeros(I,n, 2))
  end

  Base.@kwdef struct FBD{l, n, F<:AbstractFloat, I<:Integer}        
    L::I = l
    N::I = n
    x::SVector{n,F} = SVector{n}(zeros(F, n))
    neighbours::SMatrix{n, 4, I} = SMatrix{n, 4}(zeros(I, n, 4))
    neighbourhoods::SMatrix{n, 8, I} = SMatrix{n, 8}(zeros(I, n, 8))
    neighbourhood_values::SVector{n,I} = SVector{n}(zeros(I, n))

    # These values are reset for each step
    σ::SVector{n,F} = SVector{n}(zeros(F, n)) # Relative tension
    tension::SVector{n,F} = SVector{n}(zeros(F, n))
    max_σ::F = 0.0
    status::SVector{n,I} = SVector{n}(zeros(I, n))
    cluster_size::SVector{n,I} = SVector{n}(zeros(I, n))
    cluster_dimensions::SVector{n,I} = SVector{n}(zeros(I, n))
    # Relative possition of every fiber with respect to it's cluster
    rel_pos_x::SVector{n,I} = SVector{n}(zeros(I, n))
    rel_pos_y::SVector{n,I} = SVector{n}(zeros(I, n))
    cluster_outline_length::SVector{n,I} = SVector{n}(zeros(I, n))
    # These values are reset for each cluster
    cluster_outline::SVector{n,I} = SVector{n}(zeros(I, n))
    unexplored::SVector{n,I} = SVector{n}(zeros(I, n))

    # These arrays store one value for each step
    most_stressed_fiber::SVector{n,F} = SVector{n}(zeros(F, n))
    nr_clusters::SVector{n,I} = SVector{n}(zeros(I, n))
    largest_cluster::SVector{n,I} = SVector{n}(zeros(I, n))
    largest_perimiter::SVector{n,I} = SVector{n}(zeros(I, n))
end

x = Foo3{2,4, 3, Float64, Int64}()