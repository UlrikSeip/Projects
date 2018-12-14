include("RungeKutta4.jl")
include("SIR.jl")
import PyPlot
const plt = PyPlot

function main(A,a0,om,b,c,d,di,e,f_funk,fvar,S0,I0,R0,T,dt,filename)
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
        s = RungeKutta4(dt, dSall, S[j], t[j], [c,R[j],A,a0,om,I[j],N[j],d,e,f_funk,fvar])
        i = RungeKutta4(dt, dIall, I[j], t[j], [A,a0,om,S[j],N[j],b,d,di])
        r = RungeKutta4(dt, dRall, R[j], t[j], [c,I[j],b,d,S[j],f_funk,fvar])
        push!(S,s)
        push!(I,i)
        push!(R,r)
        push!(N,s+i+r)
    end
    #rs = (b/c)*(1-(b/a))/(1 + (b/c))
    print("Final tally: N=")
    println(S[end]+I[end]+R[end])
    print("Max N=")
    println(maximum(N))
    #print("Expected: ")
    #print("S=" * string(b/a) * ", I="*string((-(b/a)+1)/(1 + (b/c))))
    #println(", R=" * string(rs))
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
    #plt.savefig("C:\\Users\\Bendik\\Documents\\GitHub\\Projects\\prosjekt5\\plots/"*filename)
    #plt.savefig("/plots/"*filename)
    plt.show()
end


a0 = 6              #lingerling chem trails
A = 3               #set A = 0 to disable oscillator (size of chem trails)
om = 2*pi/356.25 #frequency of chem-trails
b = 1               #effectiveness of healing crystals, pop A
c = 0.5             #levels of atheism
d =  0.00002242299  #natural death rate
di = 0.000          #disease death rate
e =  0.00002948891  #birth rate
f0 = 0.00           #f = 0 for effective anti-vac campaigns
Df = 0.001
fvar = [f0,Df]
S0 = 300            
I0 = 100
R0 = 0
T = 365.25*5
dt = 0.01

#main(A,a0,om,b,c,d,di,e,df,fvar,S0,I0,R0,25,dt,"RK4_comp.pdf")
#main(A,a0,om,b,c,d,di,e,df,fvar,S0,I0,R0,100,dt,"RK4_comp.pdf")
main(A,a0,om,b,c,d,di,e,df,fvar,S0,I0,R0,T,dt,"RK4_comp.pdf")

