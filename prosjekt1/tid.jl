using Dates
using Statistics
include("test_b.jl")
include("test_c.jl")

n = [2, 3, 4, 5, 6, 7, 8, 9]
N = [10, Int64(1e2), Int64(1e3), Int64(1e4), Int64(1e5), Int64(1e6), Int64(1e7)]
T0 = zeros(length(N))
T1 = zeros(length(N))
m = 100

for i in range(1, stop=length(n))
    t0 = @elapsed rref_b(n[i])
    t1 = @elapsed rref_c(n[i])
end

temp = zeros(m)

for i in range(1, stop=length(N))
    for l in range(1, stop=m)
        temp[l] = @elapsed rref_b(N[i])
    end
    t0 = mean(temp)
    T0[i] = t0

    for l in range(1, stop=m)
        temp[l] = @elapsed rref_c(N[i])
    end
    t1 = mean(temp)
    T1[i] = t1
end

print("b = ")
println(T0)
print("c = ")
println(T1)

