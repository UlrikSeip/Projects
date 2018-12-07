include("RungeKutta4.jl")
include("SIR.jl")
import PyPlot
const plt = PyPlot

function ODE_sol(A,a0,om,b,c,S0,I0,R0,T,dt,filename)
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

    t = 0:dt:T+dt #an array of timesteps
    N0 = S0+I0+R0 #the total population
    N = [Float64(N0)]

    for j = 1:length(t)-1
        #finds a
        a = a_os(A,a0,om,t[j])
        #finds the new values for S,I and R
        s = RungeKutta4(dt, dS, S[j], t[j], [c,R[j],a,I[j],N[j]])
        i = RungeKutta4(dt, dI, I[j], t[j], [a,S[j],N[j],b])
        r = RungeKutta4(dt, dR, R[j], t[j], [c,I[j],b])

        push!(N,s+i+r)
        push!(S,s)
        push!(I,i)
        push!(R,r)
    end

    
    #rs = (b/c)*(1-(b/a))/(1 + (b/c))
    print("Final tally: N=")
    println(S[end]+I[end]+R[end])
    #print("Expected: ")
    #print("S=" * string(b/a) * ", I="*string((-(b/a)+1)/(1 + (b/c))))
    #println(", R=" * string(rs))
    print("Ana: S=" * string(S[end]/N[end]) * ", I=" * string(I[end]/N[end]))
    println(", R=" * string(R[end]/N[end])* ".")
    print("Number: S=" * string(S[end]) * ", I=" * string(I[end]))
    println(", R=" * string(R[end])* ".")
    println(maximum(N))
    println()
    #println("\n") 
    plt.plot(t,S)
    plt.plot(t,I)
    plt.plot(t,R)
    #plt.plot(t,N)
    plt.grid()
    plt.legend(["Susceptible", "Infected", "Recovered", "Total"])
    plt.xlabel("Days")
    plt.ylabel("People")
    #plt.savefig("C:\\Users\\Bendik\\Documents\\GitHub\\Projects\\prosjekt5\\plots/"*filename)
    #plt.savefig("/plots/"*filename)
    plt.show()
end




function a_os(A,a0,om,t)
    #we assume we start with average transmission rate
    #return A*sin.(t*om) + a0
    return A*cos.(t*om) + a0
end

T = 365.25*2
dt = 0.01
t = range(0, stop=T, step=dt)
a = zeros(length(t))

for i = 1:length(a)
    a[i] = a_os(2,4,2*pi/365.25,t[i])
end


#ODE_sol(1,4,2*pi/365.25, 1,0.5, 100,100,200, T,dt, "opp_d_A.pdf")
#ODE_sol(1,4,2*pi/365.25, 4,0.5, 100,100,200, T,dt, "opp_d_B.pdf")
#ODE_sol(2,4,2*pi/365.25, 2,0.5, 100,100,200, T,dt, "opp_d_C.pdf")
ODE_sol(2,2,2*pi/365.25, 1,0.5, 100,100,200, T,dt, "opp_d_C.pdf")

#ODE_sol(4,2,0.5, 300,100,0, 41,.01, "opp_d_B.pdf")
#ODE_sol(4,3,0.5, 300,100,0, 27,.01, "opp_d_C.pdf")
#ODE_sol(4,4,0.5, 300,100,0, 27,.01, "opp_d_D.pdf")

