include("rotator.jl")

#function for making the matrix
function make_mat(rho_min, rho_max, n, om)
    h = (rho_max-rho_min)/n #the number of steps
    e = -1/h^2 #the non-diagonal matrix element
    println(h)
    println(e)

    rho = range(rho_min, step=h, stop=rho_max) #the values of ρ
    d = zeros(n) #array to hold the diagonal matrix element
    #finds the diagonal matrix element for each ρ
    for r in range(1, stop=length(rho)-1, step=1)
        d[r] = 2/h^2 + rho[r]^2*om + 1/rho[r] #tror dette er rett
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
rho_max = 10
n = 100 #the step length
a = make_mat(rho_min, rho_max, n, 0.01)
#println(a)
new_a, r, n, counter, tol = rotate(a, 1e-8)
#god_print(r)
dia_print(new_a)
println(tol)
