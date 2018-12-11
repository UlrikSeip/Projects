include("RungeKutta4.jl")
include("SIR.jl")
import PyPlot
const plt = PyPlot

function main(a,b,c,f_funk,fvar,S0,I0,R0,T,dt,filename)
    """
    A function that runs a complete SIRS model using RK4 when accounting
    for vaccination, represented with f.
    Takes the inputs rate of transmission a, rate of recovery b, rate of
    immunity loss c, function for f f_funk, variables for f_funk fvar, 
    initial values of susceptible, infected and recovered individuals, S0, 
    I0 and R0 respectively, the total simulation time T in days, the 
    timestep dt, and file filename.
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
        s = RungeKutta4(dt, dSf, S[j], t[j], [c,R[j],a,I[j],N[j],f_funk,fvar])
        i = RungeKutta4(dt, dI, I[j], t[j], [a,S[j],N[j],b])
        r = RungeKutta4(dt, dRf, R[j], t[j], [c,I[j],b,S[j],f_funk,fvar])
        push!(S,s)
        push!(I,i)
        push!(R,r)
        push!(N,s+i+r)
    end
    rs = (b/c)*(1-(b/a))/(1 + (b/c))
    print("Final tally: N=")
    println(S[end]+I[end]+R[end])
    print("Max N=")
    println(maximum(N))
    print("Expected: ")
    print("S=" * string(b/a) * ", I="*string((-(b/a)+1)/(1 + (b/c))))
    println(", R=" * string(rs))
    print("Ana: S=" * string(S[end]/N[end]) * ", I=" * string(I[end]/N[end]))
    println(", R=" * string(R[end]/N[end])* ".")
    print("Number: S=" * string(S[end]) * ", I=" * string(I[end]))
    println(", R=" * string(R[end])* ".\n")
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

S0 = 300
I0 = 100
R0 = 0

#=
#constant f
main(4,1,0.5, df,[1,0], 300,100,0, 15,.01, "opp_e_A0.pdf")
main(4,2,0.5, df,[1,0], 300,100,0, 15,.01, "opp_e_B0.pdf")
main(4,3,0.5, df,[1,0], 300,100,0, 15,.01, "opp_e_C0.pdf")
main(4,4,0.5, df,[1,0], 300,100,0, 15,.01, "opp_e_D0.pdf")

#linear f
main(4,1,0.5, df,[0,.01], 300,100,0, 30,.01, "opp_e_A1.pdf")
main(4,1,0.5, df,[0,.01], 300,100,0, 300,.01, "opp_e_A2.pdf")
main(4,1,0.5, df,[0,.01], 300,100,0, 1500,.01, "opp_e_A3.pdf")
main(4,2,0.5, df,[0,.01], 300,100,0, 25,.01, "opp_e_B1.pdf")
main(4,2,0.5, df,[0,.01], 300,100,0, 80,.01, "opp_e_B2.pdf")
main(4,3,0.5, df,[0,.01], 300,100,0, 35,.01, "opp_e_C1.pdf")
main(4,3,0.5, df,[0,.01], 300,100,0, 60,.01, "opp_e_C2.pdf")
main(4,4,0.5, df,[0,.01], 300,100,0, 25,.01, "opp_e_D1.pdf")
main(4,4,0.5, df,[0,.01], 300,100,0, 60,.01, "opp_e_D2.pdf")
=#
#campaign f
main(4,1,0.5, f_campaign,[2,9,0,0.3], 300,100,0, 20,.01, "opp_e_fa.pdf")
main(4,2,0.5, f_campaign,[2,9,0,0.3], 300,100,0, 20,.01, "opp_e_fb.pdf")
main(4,3,0.5, f_campaign,[2,9,0,0.3], 300,100,0, 20,.01, "opp_e_fc.pdf")
main(4,4,0.5, f_campaign,[2,9,0,0.3], 300,100,0, 20,.01, "opp_e_fd.pdf")
