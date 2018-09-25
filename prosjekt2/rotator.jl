using LinearAlgebra

n=Int64(1e2)
a = ones(Float64, n, n)        #getting an error here, y?
tol = 0.0001
r = Matrix{Float64}(I, n, n)     #initialising eigenvector matrix

function god_print(noe)
    for i in range(1, step=1, length=n)
        println(noe[i,:])
    end
end

function off(a)
    return maximum(a^2) #this might work? p. 217 lecture notes.
end

function maxKnotL(a) #using this instead of ^^
    max = 0
    kl = [1,1]
    for k in range(1, step=1, length=length(a[1,:]))
        for l in range(1, step=1, length=length(a[2,:]))
            if ((a[k, l] > max) && (k != l))
                max = a[k, l]
                kl = [k, l]
            end
        end
    end
    #god_print(a)
    return kl[1], kl[2]                    #k, l                       
end

function rotate(a, r)                 #denne må doublifiseres ganske kraftig tror jeg
    b = true
    k, l = maxKnotL(a)
    while a[k ,l] > tol         
        if (a[k, l] != 0)
            kl = maxKnotL(a)
            tau = (a[l, l]-a[k, k])/2*a[k, l] #blir ikke dette alltid null?
            if (tau > 0) #kan du ikke bare skrive t = abs(1/tau) + sqrt(1+tau^2)?
                t = 1/tau + sqrt(1+tau^2)
            #elseif tau == 0     #mulig problemet ligger i at man får delt på null her
            #    println("tau er 0")
            #    t = -1/-tau + sqrt(1+tau^2)
            else
                t = -1/-tau + sqrt(1+tau^2)
            c = 1/sqrt(1+t*t)
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
        for i in range(1, step=1, length=length(a[1,:]))
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



thing, r_ = rotate(a, r)
#println(thing)
god_print(thing)


