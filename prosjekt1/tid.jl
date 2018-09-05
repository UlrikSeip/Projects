using Dates
include("test_b.jl")
include("test_c.jl")

N = [10, Int64(1e2), Int64(1e3), Int64(1e4), Int64(1e5), Int64(1e6), Int64(1e7)]
for n in N
    print("b med " * string(n) * " n: ")
    @time rref_b(n)

    print("c med " * string(n) * " n: ")
    @time rref_c(n)

end
