using Dates
using Statistics
using LinearAlgebra
include("test_c.jl")

N = [10, 100, 1000]

function lu_dec(n)
    M = zeros(n,n) #lager en matrise
    h = 1/(n+1)
    x = range(0, stop=1, length=n) #x-verdiene
    B = zeros(n) #høyre side av likningen
    B_ = zeros(n) #den reduserte B
    k = zeros(n) #det endelige resultatet

    #setter inn verdier
    M[1,1] = 2; M[1,2] = -1 #første linje
    for i in range(2, stop=n-1) #andre til nest siste linje
        M[i,i-1] = -1
        M[i,i] = 2
        M[i,i+1] = -1
        B[i] = h^2*f(x[i])
        #M[i,n+1] = h^2*f(x[i])
    end
    M[n,n-1] = -1; M[n,n] = 2 #siste linje
    lum = lu(M) #gjør lu-dekomposisjon
    L = lum.L
    U = lum.U
    B_[1] = B[1]
    #tror vi kan anta at diagonalen til L alltid er 1, satser på det
    #forward sub. av L gir B_
    for i in range(2, stop=n)
        B_[i] = B[i]-B_[i-1]*L[i,i-1]
    end
    #backward sub. av U gir svaret
    k[n] = B_[n]
    for i in range(n-1, stop=1, step=-1) 
        k[i] = B_[i] - U[i,i+1]*k[i+1]/U[i,i]
    end
    return x, k

end

function til_fil()
    res = """"""
    f = open("res_f.txt","w")
    for n in N
        x, k = lu_dec(n)
        res *= "x = " * string(x) * "\nk = " * string(k) * "\n"
    end
    #println(res)
    write(f, res)
    close(f)
end

function time_e(n)
    t = @elasped lu_dec(n)

til_fil()

