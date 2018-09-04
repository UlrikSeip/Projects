using Dates
include("test_b.jl")
include("test_c.jl")

N = [10, Int64(1e2), Int64(1e3), Int64(1e4), Int64(1e5), Int64(1e6), Int64(1e7)]
for n in N
    b_start = now()
    rref_b(n)
    b_stopp = now()

    c_start = now()
    rref_c(n)
    c_stopp = now()

    b_tid =  Microsecond(b_stopp - b_start)
    c_tid =  Microsecond(c_stopp - c_start)

    println("For n=" * string(n) * " bruker b-funksjonen " * string(b_tid) * " og c-funksjonen " * string(c_tid) * ".")
end
