include("rotator.jl")
using LinearAlgebra

#vars
n = Int64(5)
a = zeros(Float64, n, n)

for i in range(1, step = 1, length = n)
    for j in range(1, step = 1, length = n)
        if ((j == i+1) || (j == i-1))
            a[i, j] = 1                     # b = 2
        end
        if (i == j)
            a[i, j] = 2                     # a = 5
        end
    end
end

#a__ = eigvals(a)
#a_ = jacobi_method(a)
#a_, r_, n, counter, tol = rotate(a, 1e-4)
#god_print(a_)
#println(a__)
#println(find_eigvals(a_))
#aprinter()
#dimprinter()
#counterprinter()
filemaker(10, 10, 300, 1e-4, "rotated.txt")