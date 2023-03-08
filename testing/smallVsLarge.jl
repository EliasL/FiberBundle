L = 20240
N = L*L
x = Float32(1/N)
s = sum([x for _ in 1:N])
X = Float32(N)
S = sum([X for _ in 1:N])

println(s*N)
println(S/N)