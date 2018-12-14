include("RungeKutta4.jl")
include("SIR.jl")
import PyPlot
const plt = PyPlot

function ODE_sol(a,b,c,d,di,e,S0,I0,R0,T,dt,filename)
    """
    A function that runs a complete SIRS model using RK4 when accounting for
    vital dynamics, ie. people dying or being born, or dying from the 
    disease.
    Takes the inputs rate of transmission a, rate of recovery b, rate of
    immunity loss c, rate of death d, rate of death from the disease di, 
    rate of people being born e, initial values of susceptible, infected and 
    recovered individuals, S0, I0 and R0 respectively, the total simulation
    time T in days, the timestep dt, and file filename.
    Plots the resulting arrays and saves the file as filename.
    """

    S = [Float64(S0)] #lists to hold the resulting values for individuals
    I = [Float64(I0)]
    R = [Float64(R0)]
    D = [0.0] #a list to hold how much the population has decreased
    #if d=e=0 then it will display how many people have been killed by 
    #the disease

    t = 0:dt:T+dt #an array of timesteps
    N0 = S0+I0+R0 #the total population
    N = [Float64(N0)]

    for j = 1:length(t)-1
        #finds the new values for S,I and R
        s = RungeKutta4(dt, dSvd, S[j], t[j], [c,R[j],a,I[j],N[j],d,e])
        i = RungeKutta4(dt, dIvd, I[j], t[j], [a,S[j],N[j],b,d,di])
        r = RungeKutta4(dt, dRvd, R[j], t[j], [c,I[j],b,d])

        #finds the decrease in population
        Dd = D[j] -s-i-r+S[j]+I[j]+R[j] 

        push!(N,s+i+r)
        push!(D,Dd)
        push!(S,s)
        push!(I,i)
        push!(R,r)
    end

    #prints and plots the results
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
    print("Total dec, D=")
    println(D[end])
    println()
    plt.plot(t,S)
    plt.plot(t,I)
    plt.plot(t,R)
    #plt.plot(t,D)
    #plt.plot(t,N)
    plt.grid()
    plt.legend(["Susceptible", "Infected", "Recovered", "Total"]) #"Dead",
    plt.xlabel("Days")
    plt.ylabel("People")
    plt.savefig("C:\\Users\\Bendik\\Documents\\GitHub\\Projects\\prosjekt5\\plots/"*filename)
    #plt.savefig("/plots/"*filename)
    plt.show()
end

#realistic values for d and e based on Norway
d =  0.00002242299
di = 0
e =  0.00002948891

#=
#generic models with d,e and di
ODE_sol(4,1,0.5,d,di,e, 300,100,0, 15,.01, "opp_c_A0.pdf")
ODE_sol(4,1,0.5,d,0.1,e, 300,100,0, 25,.01, "opp_c_A1.pdf")
ODE_sol(4,2,0.5,d,di,e, 300,100,0, 25,.01, "opp_c_B0.pdf")
ODE_sol(4,2,0.5,d,0.1,e, 300,100,0, 30,.01, "opp_c_B1.pdf")
ODE_sol(4,3,0.5,d,di,e, 300,100,0, 35,.01, "opp_c_C0.pdf")
ODE_sol(4,3,0.5,d,.1,e, 300,100,0, 35,.01, "opp_c_C1.pdf")
ODE_sol(4,4,0.5,d,di,e, 300,100,0, 27,.01, "opp_c_D0.pdf")
ODE_sol(4,4,0.5,d,.1,e, 300,100,0, 27,.01, "opp_c_D1.pdf")
=#
#different attempts to kill the entire population
#ODE_sol(4,1,0.5, d,3,e, 300,100,0, 60,.01, "opp_c_k0.pdf")
#ODE_sol(4,1,0.5, d,5,e, 300,100,0, 50,.01, "opp_c_k1.pdf")
#ODE_sol(4,1,0.5, d,.1,e, 300,100,0, 50,.01, "opp_c_k2.pdf")
#ODE_sol(4,1,0.5, d,.1,e, 300,100,0, 350,.01, "opp_c_k2l.pdf")
ODE_sol(4,1,0.5, d,.01,e, 300,100,0, 350,.01, "opp_c_k3.pdf")
#=
#the growth when d and e are larger
mod = 1e3
ODE_sol(4,2,0.5, d*mod,.05,e*mod, 300,100,0, 50,.01, "opp_c_h_"*string(mod)*".pdf")
mod = 2e3
ODE_sol(4,2,0.5, d*mod,.05,e*mod, 300,100,0, 50,.01, "opp_c_h_"*string(mod)*".pdf")
mod = 1e4
ODE_sol(4,2,0.5, d*mod,.05,e*mod, 300,100,0, 50,.01, "opp_c_h_"*string(mod)*".pdf")
=#
#ODE_sol(4,2,0.5,d,di,e, 300,100,0, 365.25*10,.01, "opp_c_B0.pdf")
