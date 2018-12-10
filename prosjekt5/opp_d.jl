include("RungeKutta4.jl")
include("SIR.jl")
import PyPlot
const plt = PyPlot

function ODE_sol(A,a0,om,b,c,S0,I0,R0,T,dt,filename)
    """
    A function that runs a complete SIRS model using RK4, when the rate of 
    transmission a can vary with time. 
    Takes the inputs average transmission rate a0, maximum deviation from 
    a0 A, the frequency of oscillation om, rate of recovery b, rate of
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
        #finds the new values for S,I and R
        s = RungeKutta4(dt, dSda, S[j], t[j], [c,R[j],A,a0,om,I[j],N[j]])
        i = RungeKutta4(dt, dIda, I[j], t[j], [A,a0,om,S[j],N[j],b])
        r = RungeKutta4(dt, dR, R[j], t[j], [c,I[j],b])

        push!(N,s+i+r)
        push!(S,s)
        push!(I,i)
        push!(R,r)
    end

    #prints and plots the results
    print("Final tally: N=")
    println(S[end]+I[end]+R[end])
    print("Max N=")
    println(maximum(N))
    print("Ana: S=" * string(S[end]/N[end]) * ", I=" * string(I[end]/N[end]))
    println(", R=" * string(R[end]/N[end])* ".")
    print("Number: S=" * string(S[end]) * ", I=" * string(I[end]))
    println(", R=" * string(R[end])* ".")
    println()
    plt.plot(t,S)
    plt.plot(t,I)
    plt.plot(t,R)
    #plt.plot(t,N)
    plt.grid()
    plt.legend(["Susceptible", "Infected", "Recovered", "Total"])
    plt.xlabel("Days")
    plt.ylabel("People")
    plt.savefig("C:\\Users\\Bendik\\Documents\\GitHub\\Projects\\prosjekt5\\plots/"*filename)
    #plt.savefig("/plots/"*filename)
    plt.show()
end

T = 365.25*2 
dt = 0.01

#if om=2*pi/365.25 a will have cycle of one year
ODE_sol(1,4,2*pi/365.25, 1,0.5, 100,100,200, T,dt, "opp_d_A.pdf")
ODE_sol(1,4,2*pi/365.25, 4,0.5, 100,100,200, T,dt, "opp_d_B.pdf")
ODE_sol(2,4,2*pi/365.25, 2,0.5, 100,100,200, T,dt, "opp_d_C.pdf")
ODE_sol(2,6,2*pi/365.25, 2,0.5, 100,100,200, T,dt, "opp_d_D.pdf")


