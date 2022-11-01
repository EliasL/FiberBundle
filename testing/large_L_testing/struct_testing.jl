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


x = Foo3{2,4, 3, Float64, Int64}()
x = A(c = 4)