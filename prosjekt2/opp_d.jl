include("rotator.jl")

#function for making the matrix
function make_mat(rho_min, rho_max, n)
    h = (rho_max-rho_min)/n #the number of steps
    e = -1/h^2 #the non-diagonal matrix element
    println(h)
    println(e)

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
    #god_print(a)
    return a
end

rho_min = 0
<<<<<<< HEAD
rho_max = 5
h = 1. #the step length
a = make_mat(rho_min, rho_max, h)
=======
rho_max = 10
n = 100 #the step length
a = make_mat(rho_min, rho_max, n)
>>>>>>> 209008b159fa7d5a46e05d8a4b32a5acbf6c5dca
#println(a)
god_print(a)

new_a, r, n, counter, tol = rotate(a, 1e-8)
<<<<<<< HEAD
#god_print(new_a)
god_print(a)
=======
#god_print(r)
>>>>>>> 209008b159fa7d5a46e05d8a4b32a5acbf6c5dca
dia_print(new_a)
println(tol)
