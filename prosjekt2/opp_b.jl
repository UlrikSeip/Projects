include("rotator.jl")

#vars
n = Int64(4)
a = ones(Float64, n, n)

#a_, r_, n, counter, tol = rotate(a, 1e-10)
#eigenprinter()
#aprinter()
#dimprinter()
#counterprinter()
filemaker(5, 5, 200, 1e-10, "rotated.txt")