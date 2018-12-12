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

function MC_sol(a,b,c,S0,I0,R0,T,sims,filename)
    """
    Takes the inputs rate of transmission a, rate of recovery b, rate of
    immunity loss c, initial values of susceptible, infected and 
    recovered individuals, S0, I0 and R0 respectively, the total simulation
    time T in days, the timestep dt, and file filename.
    Plots the resulting arrays and saves the file as filename.
    """

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
    @time Threads.@threads for j = 1:sims

        N = [Float64(N0)]
        S = [Float64(S0)] #lists to hold the resulting values for individuals
        I = [Float64(I0)]
        R = [Float64(R0)]

        for i = 1:length(t)-1
            #finds the new values for S,I and R
            s, i, r = mcStep(a, b, c, S[i], I[i], R[i], N[i], dt)
            #finds the decrease in population
            #Dd = D[j] -s-i-r+S[j]+I[j]+R[j]
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
    Navg, Nstd = statisticsinator(Ns)
    Savg, Sstd = statisticsinator(Ss)
    Iavg, Istd = statisticsinator(Is)
    Ravg, Rstd = statisticsinator(Rs)

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
    println("Generating plots")
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

MC_sol(4, 1, 0.5, 300, 100, 0, 20, 10000, "b.pdf")