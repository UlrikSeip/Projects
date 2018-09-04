include("test_c.jl")

N = [10, Int64(1e2), Int64(1e3), Int64(1e4), Int64(1e5), Int64(1e6), Int64(1e7), Int64(1e8)]

#går gjennom alle n-ene
for n in N
    x, k, b_, u = rref_c(n) #kaller funksjonen

    eps_max = abs((k[n]-u[n]/u[n])) #maks feil
    #finner feilen i alle punktene, hvis feilen er større enn maks feil blir den maks feil
    for i in range(1, stop=length(k)-1, step=1)
        eps = abs((k[i]-u[i]/u[i]))
        if eps > eps_max
            eps_max = eps
        end
    end
    eps_max = abs(log10(eps_max))
    print(string(eps_max) * ", ")
end
