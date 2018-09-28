include("rotator.jl")

#function for making the matrix
function make_mat(rho_min, rho_max, h)
    n = Int64((rho_max-rho_min)/h) #the number of steps
    e = -1/h^2 #the non-diagonal matrix element

    rho = range(rho_min, step=h, stop=rho_max) #the values of ρ
    d = zeros(n) #array to hold the diagonal matrix element
    #finds the diagonal matrix element for each ρ
    for r in range(1, stop=length(rho)-1, step=1)
        d[r] = 2/h^2 + rho[r]^2
    end

    a = zeros(Float64, n, n) #the matrix, we now have to fill in the diagonal values
    a[1,1] = d[1] #the first elements
    a[1,1+1] = e
    for i in range(2, stop=n-1, step=1) #all the other elements
        a[i,i-1] = e
        a[i,i] = d[i]
        a[i,i+1] = e
    end
    a[n,n-1] = e #the last elements
    a[n,n] = d[n]
    return a
end

rho_min = 0
rho_max = 1e2
h = 1. #the step length
a = make_mat(rho_min, rho_max, h)
#println(a)
new_a, r, n, counter, tol = rotate(a, 1e-8)
#god_print(new_a)
god_print(r)
println(tol)

