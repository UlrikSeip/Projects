include("rotator.jl")
include("project2.jl")

#function for making the matrix
function make_mat(rho_minn, rho_max, n)
    h = (rho_max-rho_minn)/n #the number of steps
    e = -1/h^2 #the non-diagonal matrix element

    d = zeros(n) #array to hold the diagonal matrix element
    #finds the diagonal matrix element for each ρ
    #d[1] = 2/h^2 + rho_minn^2
    for r=1:n
        d[r] = 2/h^2 + (rho_minn + r*h)^2
    end

    a = zeros(Float64, n, n) #the matrix, we now have to fill in the diagonal values
    a[1,1] = 2/h^2 + rho_minn^2 #the first elements
    a[1,1+1] = e
    for i=2:n-1 #all the other elements
        a[i,i-1] = e
        a[i,i] = d[i]
        a[i,i+1] = e
    end
    a[n,n-1] = e #the last elements
    a[n,n] = 2/h^2 + rho_max^2
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

#a function that finds the eigenvalues analyticaly
function ana_eigen(n)
    eig = zeros(n)
    #some constants
    hbar = 4.135667662e-15 
    omega = 1
    m = 0.5110e-6
    k = m*omega^2
    alpha = (hbar^2/(m*k))^0.25
    l_factor = (2*m*alpha^2)/(hbar^2)

    for i=1:n
        eig[i] = l_factor*hbar*omega*(2*(i-1) + 1.5)
    end
    return sort(eig)
end

#the main part of the task, we go from rho_min to rho_max with n itterations
function main(rho_min, rho_max, n)
    a = make_mat(rho_min, rho_max, n)
    #a = SQM_Tridiag(n)

    new_a, r, n, counter, tol = rotate(a, 1e-4)
    dia = find_dia(new_a)
    dia = sort(dia)
    print(string(rho_max) * ": ")

    eigen = ana_eigen(n)
    diff = dia-eigen
    l_mean = stat.mean(diff)
    l_std = stat.std(diff)
    println("l_mean: " * string(l_mean) * " l_std: " * string(l_std))
    println(dia[1:20])
    #println(eigen[1:20])
    println()
    return dia, eigen
end



"""
main(1e-6, 1, 100)
main(1e-6, 10, 100)
main(1e-6, 20, 100)
main(1e-6, 30, 100)
main(1e-6, 40, 100)
main(1e-6, 50, 100)
main(1e-6, 60, 100)
main(1e-6, 70, 100)
main(1e-6, 80, 100)
main(1e-6, 90, 100)
main(1e-6, 100, 100)


main(1e-6, 11, 100)
main(1e-6, 12, 100)
main(1e-6, 13, 100)
main(1e-6, 14, 100)
main(1e-6, 15, 100)
main(1e-6, 16, 100)
main(1e-6, 17, 100)
main(1e-6, 18, 100)
main(1e-6, 19, 100)
main(1e-6, 21, 100)
main(1e-6, 22, 100)
main(1e-6, 23, 100)
main(1e-6, 24, 100)
main(1e-6, 25, 100)
main(1e-6, 26, 100)
main(1e-6, 27, 100)
main(1e-6, 28, 100)
main(1e-6, 29, 100)
"""

#main(1e-6, 11, 200)
#main(1e-6, 11, 300)
#main(1e-6, 11, 400)
#main(1e-6, 10, 400)
#main(1e-6, 11, 500)
#main(1e-6, 9, 400)

"""
main(1e-6, 8, 400)
main(1e-6, 10, 400)
main(1e-6, 12, 400)
main(1e-6, 13, 400)
main(1e-6, 14, 400)
"""


main(1e-6, 4, 400)
main(1e-6, 5, 400)
#main(1e-6, 6, 400)
#main(1e-6, 7, 400)
