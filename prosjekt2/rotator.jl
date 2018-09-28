using LinearAlgebra

#vars
"""
n = Int64(20)
a = ones(Float64, n, n)
"""
tol = 1e-10


"""
call rotate(a) to run
outputs newA, r, n, counter
print outputs can be de-commented on the final lines of the script
"""
#prints a matrix in a somewhat readable format
function god_print(noe)
    for i in range(1, step=1, length=n)
        println(noe[i,:])
    end
end

#prints the diagonal of a matrix "noe"
function dia_print(noe)
    for i in range(1, step=1, length=n)
        print(noe[i,i])
        print(", ")
    end
    println(" ")
end

function off(a)
    return maximum(a^2) #this might work? p. 217 lecture notes.
end

function maxKnotL(a) #using this instead of ^^
    max = 0
    kl = [1,1]
    n = Int64(length(a[1,:]))
    for k in range(1, step=1, length = n)
        for l in range(1, step=1, length = n)
            if ((a[k, l] > max) && (k != l))
                max = a[k, l]
                kl = [k, l]
            end
        end
    end
    #god_print(a)
    return kl[1], kl[2]                    #k, l                       
end

function rotate(a)  
    n = Int64(length(a[1,:]))        #initiates n for later use
    r = Matrix{Float64}(I, n, n)     #initialising eigenvector matrix
    counter = 0                      #teller antall "similarity transformaitons"
    k, l = maxKnotL(a)               #finds indices of matrix element with highest value  
    while a[k ,l] > tol              #this is the actual loop
        counter += 1                 #
        if (a[k, l] != 0)
            kl = maxKnotL(a)
            tau = (a[l, l]-a[k, k])/2*a[k, l] #blir ikke dette alltid null?
            #jupp, fuck all this
            """
            #if (tau > 0) #kan du ikke bare skrive t = abs(1/tau) + sqrt(1+tau^2)?
                #t = 1/tau + sqrt(1+tau^2)
            #    t = tau + sqrt(1+tau^2)
            #elseif tau == 0     #mulig problemet ligger i at man får delt på null her
            #    println("tau er 0")
            #    t = -1/-tau + sqrt(1+tau^2)
            #else
                #t = -1/-tau + sqrt(1+tau^2)
            #    t = tau + sqrt(1+tau^2)
            """
            t = -tau + sqrt(1+tau^2)
            c = 1/sqrt(1+t^2)
            s = c*t
        #end
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
    #god_print(r)

    return a, r, n, counter
end

#formatted printing functions:

function eigenprinter()
    println("Egienvectors:")
    println()
    god_print(r_)
    println()
end

function aprinter()
    println("Transformed matrix:")
    println()
    god_print(a_)
    println()
end

function dimprinter()
    println("Matrix dim: ", n, "*", n)
    println()
end

function counterprinter()
    println("Number of similarity transformations: ", counter)
end

function filemaker(start, stop)                             #creates rotated.txt with format "n, counter \n"
    open("rotated.txt", "w") do f                           #clears file
    end
    for i in range(start, step = 10, length = stop)         #does the actual writing from n = start to n = stop
        n = i
        a = ones(Float64, n, n)
        a_, r_, n, counter = rotate(a)
        open("rotated.txt", "a") do f
            write(f,string(n, " ", counter, "\n"))
        end
    end
end

#a_, r_, n, counter = rotate(a)
#eigenprinter()
#aprinter()
#dimprinter()
#counterprinter()

filemaker(10, 1000)