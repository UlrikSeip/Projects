#en funksjon for f(x)
function f(x)
    return (100*exp(-10*x))
end

#funksjon for radreduksjon og ting
function rref_c(n) 
    a = -1 #linjen under diagonalen
    b = 2 #diagonalen
    c = -1 #linjen over diagonalen

    B = zeros(n) #høyre side av likningen
    b_ = zeros(n) #den reduserte b
    B_ = zeros(n) #den reduserte B
    k = zeros(n) #det endelige resultatet
    x = range(0, stop=1, length=n) #x-verdiene
    
    l = zeros(n)
    for i in range(1, length=n, step=1)
        l[i] = exp(-10.0*x[i])
    end
    u = fill(1.0, n) - (1.0 - exp(-10.0))*x - l #den eksakte løsningen

    h = 1/(n+1) #steglengde
    for i in range(1, length=n, step=1) # setter inn verdiene til B
        B[i] = h^2*f(x[i])
    end 

    b_[1] = b
    B_[1] = B[1]

    for i in range(2, stop=n, step=1) #forward substitution, 4 FLOPS
        b_[i] = 2 - 1/b_[i-1]
        B_[i] = B[i] + B_[i-1]/b_[i-1]
    end

    k[n] = B_[n]
    for i in range(n-1, stop=1, step=-1) #backward substitution, 3 FLOPS
        k[i] = B_[i] + k[i+1]/b_[i] #-*(-1)
    end
    return x, k, b_, u

end
