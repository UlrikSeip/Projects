using LinearAlgebra

n=Int(1e2)
a = ones(Float64, n, n)        #getting an error here, y?
tol = 1e-8
r = Matrix{Float64}(I, n, n)     #initialising eigenvector matrix

function off(a)
    return maximum(a^2) #this might work? p. 217 lecture notes.
end

function maxKnotL(a) #using this instead of ^^
    max = 0
    for k in range(length(a[1]))
        for l in range(length(a[2]))
            if ((a[k, l] > max) && (k != l))
                max = a[k, l]
                kl = (k, l)
            end
        end
    end
    return kl[0], kl[1]                    #k, l                       
end

function rotate(a, r)                 #denne må doublifiseres ganske kraftig tror jeg
    b = true
    k, l = maxKnotL(a)
    while a[k ,l] > tol         
        if (a[k, l] != 0)
            kl = maxKnotL(a)
            tau = (a[l, l]-a[k, k])/2*a[k, l]
            if (tau > 0)
                t = 1/tau + sqrt(1+tau^2)
            else
                t = -1/-tau + sqrt(1+tau^2)
            c = 1/sqet(1+t*t)
            s = c*t
            end
        else
            c = 1
            s = 0
        end
        a_kk = a[k, k]
        a_ll = a[l, l]
        #doing stuff with indices k and l
        a[k, k] = c*c*a_kk - 2*c*s*a[k, l] + s^2*a_ll
        a[l, l] = s*s*a_kk + 2*s*s*a[k, l] + c^2*a_ll
        a[k, l] = 0
        a[l, k] = 0
        #doing stuff with the remaing matrix elements
        for i in range(Int64, len(a[0]))
            if (i != k && i != l)
                a_ik = a[i, k]
                a_il = a[i, l]
                a[i, k] = c*a_ik - s*a_il
                a[k, i] = a[i, k]
                a[i, l] = c*a_il + s*a_ik
                a[l, i] = a[i, l]
            end
               #calculating eigenvectors
            r_ik = r[i, k]
            r_il = r[i, l]
            r[i, k] = c*r_ik - s*r_il
            r[i, l] = c*r_il + s*r_ik
        end
        k, l = maxKnotL(a)
    end
    return a, r
end

thing = rotate(a, r)
print(thing)
