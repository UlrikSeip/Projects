include("rotator.jl")

#function for making the matrix
function make_mat(rho_min, rho_max, n, om)
    h = (rho_max-rho_min)/n #the number of steps
    e = -1/h^2 #the non-diagonal matrix element

    d = zeros(n) #array to hold the diagonal matrix element
    #finds the diagonal matrix element for each ρ
    for r=1:n
        d[r] = 2/h^2 + (rho_min + r*h)^2*om^2 + 1/(rho_min + r*h) #tror dette er rett
    end

    a = zeros(Float64, n, n) #the matrix, we now have to fill in the diagonal values
    a[1,1] = 2/h^2 + rho_min^2*om^2 + 1/rho_min #the first elements
    a[1,1+1] = e
    for i in range(2, stop=n-1, step=1) #all the other elements
        a[i,i-1] = e
        a[i,i] = d[i]
        a[i,i+1] = e
    end
    a[n,n-1] = e #the last elements
    a[n,n] = 2/h^2 + rho_max^2*om^2 + 1/rho_max
    return a
end

#a function theat finds the diagonal values of an matrix
function find_dia(a)
    n = length(a[1,:])
    res = zeros(n)
    for i=1:n
        res[i] = a[i,i]
    end
    return res
end

#the main part of the task, we go from rho_min to rho_max with n itterations
function main(rho_min, rho_max, n, om)
    a = make_mat(rho_min, rho_max, n, om)

    new_a, r, n, counter, tol = rotate(a, 1e-4)
    dia = find_dia(new_a)
    dia = sort(dia)
    println("rho_max " * string(rho_max) * " ohm " * string(om) * ": ")

    println(dia)
    println()
    return dia
end

main(1e-6, 23, 400, 0.01)
main(1e-6, 23, 400, 0.5)
main(1e-6, 23, 400, 1)
main(1e-6, 23, 400, 5)
