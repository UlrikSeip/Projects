include("rotator.jl")

#vars
n = Int64(5)
a = zeros(Float64, n, n)
for i in range(1, step = 1, length = n)
    for j in range(1, step = 1, length = n)
        if ((j = i+1) || j = i-1)
            a[i, j] = 2                     # b = 2
        end
        if (i = j)
            a[i, j] = 5                     # a = 5
        end
    end
end

a_, r_, n, counter, tol = rotate(a, 1e-10)
eigenprinter()
aprinter()
#dimprinter()
#counterprinter()
filemaker(5, 5, 300, 1e-10)         #takes (startN, stepN, stopN, tol)