#en funksjon for f(x)
function f(x)
    return (100*exp(-10*x))
end

#funksjon for radreduksjon og ting
function rref(n) 
    a = fill(-1, n-1) #linjen under diagonalen
    b = fill(2, n) #diagonalen
    c = fill(-1, n-1) #linjen over diagonalen

    B = zeros(n) #høyre side av likningen
    b_ = zeros(n) #den reduserte b
    B_ = zeros(n) #den reduserte B
    k = zeros(n) #det endelige resultatet
    x = range(0, stop=1, length=n) #x-verdiene
    #u = 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x) #den eksakte løsningen, fungerer ikke av en eller annen grunn

    h = 1/(n+1)
    for i in range(1, length=n, step=1)
        B[i] = h^2*f(x[i])
    end 

    b_[1] = b[1]
    B_[1] = B[1]

    for i in range(2, stop=n, step=1)
        temp = a[i-1]/b_[i-1]
        b_[i] = b[i] - temp*c[i-1]
        B_[i] = B[i] - temp*B_[i-1]
    end

    k[n] = B_[n]
    for i in range(n-1, stop=1, step=-1)
        k[i] = B_[i] - c[i]*k[i+1]/b_[i]
    end
    return x, k, b_#, u

end

rref(3)
println(2)