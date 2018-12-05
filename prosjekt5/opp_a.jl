include("RungeKutta4.jl")
include("SIR.jl")
import PyPlot
const plt = PyPlot

function main(a,b,c,S0,I0,R0,T,dt,filename)
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
    N = S0+I0+R0 #the total population

    for j = 1:length(t)-1
        #N = S[j] + I[j] + R[j]
        #= if N > 401
            println("For høy N")
            println(j)
        end
         =#
        #finds the new values for S,I and R
        s = RungeKutta4(dt, dS, S[j], t[j], [c,R[j],a,I[j],N])
        i = RungeKutta4(dt, dI, I[j], t[j], [a,S[j],N,b])
        r = RungeKutta4(dt, dR, R[j], t[j], [c,I[j],b])
        push!(S,s)
        push!(I,i)
        push!(R,r)
    end
    print("Final tally: N=")
    print(S[end]+I[end]+R[end])
    println(", S=" * string(S[end]) * ", I=" * string(I[end]) * ", R=" * string(R[end])* ".\n")
    plt.plot(t,S)
    plt.plot(t,I)
    plt.plot(t,R)
    plt.grid()
    plt.legend(["Susceptible", "Infected", "Recovered"])
    plt.xlabel("Days")
    plt.ylabel("People")
    plt.savefig("C:\\Users\\Bendik\\Documents\\GitHub\\Projects\\prosjekt5\\plots/"*filename)
    #plt.savefig("/plots/"*filename)
    plt.show()
end

main(4,1,0.5, 300,100,0, 15,.01, "opp_a_A.pdf")
main(4,2,0.5, 300,100,0, 21,.01, "opp_a_B.pdf")
main(4,3,0.5, 300,100,0, 27,.01, "opp_a_C.pdf")
main(4,4,0.5, 300,100,0, 27,.01, "opp_a_D.pdf")