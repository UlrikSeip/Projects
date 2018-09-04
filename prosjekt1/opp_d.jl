include("test_c.jl")

N = [10, Int64(1e2), Int64(1e3), Int64(1e4), Int64(1e5), Int64(1e6), Int64(1e7)]

#går gjennom alle n-ene
for n in N
    x, k, b_, u = rref(n) #kaller funksjonen

    eps_max = 0 #maks feil
    #finner feilen i alle punktene, hvis feilen er større enn maks feil blir den maks feil
    for i in range(1, stop=length(k), step=1)
        eps = log10(abs((k[i]-u[i]/u[i])))
        if eps > eps_max
            eps_max = eps
        end
    end
    println(eps_max)
end
