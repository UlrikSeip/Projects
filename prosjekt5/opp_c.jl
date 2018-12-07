include("RungeKutta4.jl")
include("SIR.jl")
import PyPlot
const plt = PyPlot

function ODE_sol(a,b,c,d,di,e,S0,I0,R0,T,dt,filename)
    """
    A function that runs a complete SIRS model using RK4.
    Takes the inputs rate of transmission a, rate of recovery b, rate of
    immunity loss c, initial values of susceptible, infected and recovered
    individuals, S0, I0 and R0 respectively, the total simulation time T in
    days, the timestep dt, and file filename.
    Plots the resulting arrays and saves the file as filename.
    """

    S = [Float64(S0)] #lists to hold the resulting values for individuals
    I = [Float64(I0)]
    R = [Float64(R0)]
    D = [0.0]
    dead = 0.0

    t = 0:dt:T+dt #an array of timesteps
    N0 = S0+I0+R0 #the total population
    N = [Float64(N0)]

    for j = 1:length(t)-1
        #N = S[j] + I[j] + R[j]
        #= if N > 401
            println("For høy N")
            println(j)
        end
         =#
        #finds the new values for S,I and R
        s = RungeKutta4(dt, dSvd, S[j], t[j], [c,R[j],a,I[j],N[j],d,e])
        i = RungeKutta4(dt, dIvd, I[j], t[j], [a,S[j],N[j],b,d,di])
        r = RungeKutta4(dt, dRvd, R[j], t[j], [c,I[j],b,d])

        Dd = D[j] -s-i-r+S[j]+I[j]+R[j] #this dosen't work if we include  
                                        #normal vital dynamics
        #dead += -s-i-r+S[j]+I[j]+R[j]
        push!(N,s+i+r)
        push!(D,Dd)
        push!(S,s)
        push!(I,i)
        push!(R,r)
    end

    rs = (b/c)*(1-(b/a))/(1 + (b/c))
    print("Final tally: N=")
    println(S[end]+I[end]+R[end])
    print("Expected: ")
    print("S=" * string(b/a) * ", I="*string((-(b/a)+1)/(1 + (b/c))))
    println(", R=" * string(rs))
    print("Ana: S=" * string(S[end]/N[end]) * ", I=" * string(I[end]/N[end]))
    println(", R=" * string(R[end]/N[end])* ".")
    print("Number: S=" * string(S[end]) * ", I=" * string(I[end]))
    println(", R=" * string(R[end])* ".")
    println(D[end])
    #println(dead)
    #println("\n")
    plt.plot(t,S)
    plt.plot(t,I)
    plt.plot(t,R)
    #plt.plot(t,D)
    #plt.plot(t,N)
    plt.grid()
    plt.legend(["Susceptible", "Infected", "Recovered", "Total"]) #"Dead",
    plt.xlabel("Days")
    plt.ylabel("People")
    #plt.savefig("C:\\Users\\Bendik\\Documents\\GitHub\\Projects\\prosjekt5\\plots/"*filename)
    #plt.savefig("/plots/"*filename)
    plt.show()
end

d =  0.00002242299
#di = 0.01
di = 0
e =  0.00002948891

#ODE_sol(4,1,0.5,d,di,e, 300,100,0, 25,.01, "opp_c_A.pdf")
#ODE_sol(4,1,0.5,d,0.1,e, 300,100,0, 230,.01, "opp_c_A.pdf")
#ODE_sol(4,2,0.5,d,di,e, 300,100,0, 41,.01, "opp_c_B.pdf")
#ODE_sol(4,2,0.5,d,0.1,e, 300,100,0, 41,.01, "opp_c_B.pdf")
#ODE_sol(4,3,0.5,d,di,e, 300,100,0, 27,.01, "opp_c_C.pdf")
#ODE_sol(4,4,0.5,d,di,e, 300,100,0, 27,.01, "opp_c_D.pdf")

#ODE_sol(4,1,0.5, 0,.1,0, 300,100,0, 350,.01, "opp_c_A.pdf")
ODE_sol(4,2,0.5, d*1,.05,e*1, 300,100,0, 35,.01, "opp_c_A.pdf")


