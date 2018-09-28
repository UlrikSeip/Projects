using LinearAlgebra


"""
call rotate(a) to run
outputs newA, r, n, counter
print outputs can be de-commented on the final lines of the script
"""
#prints a matrix in a somewhat readable format
function god_print(noe)
    for i in range(1, step=1, length = n)
        println(noe[i,:])
    end
end

function theos_print(noe)
    for i in range(1, step=1, length = n)
        println(noe[:, i])
    end
end

#prints the diagonal of a matrix "noe"
function dia_print(noe)
    for i in range(1, step=1, length = n)
        print(noe[i,i])
        print(", ")
    end
    println(" ")
end

function off(a)
    return maximum(a^2) #this might work? p. 217 lecture notes.
end

function maxKnotL(a)
    max = 0
    kl = [1,1]
    n = Int64(length(a[1,:]))
    for k in range(1, step=1, length = n)
        for l in range(1, step=1, length = n)
            if ((abs(a[k, l]) > max) && (k != l))
                max = abs(a[k, l])
                kl = [k, l]
            end
        end
    end
    return kl[1], kl[2]                    #k, l                       
end

function rotate(a, tol)              #den faktiske rotasjonsløkken. Tar inn en matrise a og nøyaktighet.
    n = Int64(length(a[1,:]))        #initiates n for later use
    r = Matrix{Float64}(I, n, n)     #initialising eigenvector matrix
    counter = 0                      #teller antall "similarity transformaitons"
    k, l = maxKnotL(a)               #finds indices of matrix element with highest value  
    while abs(a[k ,l]) > tol         #this is the actual loop
        counter += 1                 #
        if (a[k, l] != 0.0)
            #kl = maxKnotL(a)
            tau = (a[l, l] - a[k, k])/2*a[k, l] #blir ikke dette alltid null?
            if (tau > 0)
                t = -tau + sqrt(1+tau^2)
            else
                t = -tau - sqrt(1+tau^2)
            end
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
        c2 = c^2
        s2 = s^2
        a[k, k] = c2*a_kk - 2.0*c*s*a[k, l] + s2*a_ll
        a[l, l] = s2*a_kk + 2.0*c*s*a[k, l] + c2*a_ll
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
        #println(a[k, l])
    end
    #god_print(r)

    return a, r, n, counter, tol         #a = newA, r = egenvektorer, n = dim(a), counter = n(simTrans) og tol = max value for ikkediagonale elementer
end

#formatted printing functions:

function eigenprinter()
    println("Egienvectors:")
    println()
    theos_print(r_)
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

function filemaker(start, step, stop, tol, filename)        #creates filename with format "n, counter \n" and name rotated.txt
    open("rotated.txt", "w") do f                           #clears file
    end
    length = Int64(stop/step)
    for i in range(start, step = step, length = length)         #does the actual writing from n = start to n = stop
        n = i
        a = ones(Float64, n, n)
        data, time = @timed rotate(a, tol)
        a_, r_, n, counter = data[1], data[2], data[3], data[4]
        open(filename, "a") do f
            write(f,string(n, " ", counter, " ", time, "\n"))
        end
    end
end

