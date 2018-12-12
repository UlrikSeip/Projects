include("SIR.jl")
include("MC_SIR_step.jl")
include("utils.jl")
println("Importing things")
@time using PyPlot
@time import Base.Threads
@time using Statistics
const plt = PyPlot

const threads = Sys.CPU_THREADS
Threads.nthreads() = threads

function geta(t, a0, A, omega)
    return A*cos(omega*t) + a0
end

function getf(t, f0, F, omegaf)
    F*cos(omegaf*t) + f0
end

function MC_sol(a0, A, omega, b, c, d, di, bi, f0, F, omegaf, S0, I0, R0, T, sims, filename)
    """
    Takes the inputs rate of transmission a, rate of recovery b, rate of
    immunity loss c, initial values of susceptible, infected and 
    recovered individuals, S0, I0 and R0 respectively, the total simulation
    time T in days, the timestep dt, and file filename.
    Plots the resulting arrays and saves the file as filename.
    """
    a = geta(0, a0, A, omega)
    N0 = S0+I0+R0 #the total population
    tmins = [4/(a*N0), 1/(b*N0), 1/c*N0]
    dt = minimum(tmins)
    t = 0:dt:T+dt #an array of timesteps
    nts = length(t) #number of timesteps
    Ns = []
    Ss = []
    Is = []
    Rs = []

    println("Starting simulation")
    @time @inbounds Threads.@threads for j = 1:sims

        N = [Float64(N0)]
        S = [Float64(S0)] #lists to hold the resulting values for individuals
        I = [Float64(I0)]
        R = [Float64(R0)]

        @inbounds for q = 1:length(t)-1
            #new a in case of oscillation
            a = geta(t[q], a0, A, omega)
            f = getf(t[q], f0, F, omegaf)
            #finds the new values for S,I,R and N
            s, i, r = mcStep(A, omega, t[q], a0, b, c, d, di, bi, f, S[q], I[q], R[q], N[q], dt)
            push!(N,s+i+r)
            push!(S,s)
            push!(I,i)
            push!(R,r)
        end
        push!(Ns, N)
        push!(Ss, S)
        push!(Is, I)
        push!(Rs, R)
    end
    
    println("Extracting data")
    @time Navg, Nstd = statisticsinator(Ns)
    @time Savg, Sstd = statisticsinator(Ss)
    @time Iavg, Istd = statisticsinator(Is)
    @time Ravg, Rstd = statisticsinator(Rs)

    #prints and plots the results
    rs = (b/c)*(1-(b/a))/(1 + (b/c))

    print("Final tally: N=")
    println(Savg[end]+Iavg[end]+Ravg[end])
    print("Expected: ")
    print("S=" * string(b/a) * ", I="*string((-(b/a)+1)/(1 + (b/c))))
    println(", R=" * string(rs))
    print("Found: S=" * string(Savg[end]/Navg[end]) * ", I=" * string(Iavg[end]/Navg[end]))
    println(", R=" * string(Ravg[end]/Navg[end])* ".")
    print("Number: S=" * string(Savg[end]) * ", I=" * string(Iavg[end]))
    println(", R=" * string(Ravg[end])* ".")
    #print("Total dec, D=")
    #println(D[end])
    println("Standard deviations: S=$Sstd, I=$Istd, R=$Rstd")
    println()
    plt.plot(t,Savg)
    plt.plot(t,Iavg)
    plt.plot(t,Ravg)
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

a0 = 4              #lingerling chem trails
A = 3               #set A = 0 to disable oscillator (size of chem trails)
omega = 1/100       #frequency of chem-trails
b = 0.5             #effectiveness of healing crystals
c = 0.5             #levels of atheism
d = 0.00002242299
di = 0
bi = 0.00002948891
f0 = 0.001           #f = 0 for effective anti-vac campaigns
F = 1
omegaf = 1/364
S0 = 300            
I0 = 100
R0 = 0
T = 1000
sims = 1000

MC_sol(a0, A, omega, b, c, d, di, bi, f0, F, omegaf, S0, I0, R0, T, sims, "MC.pdf")
