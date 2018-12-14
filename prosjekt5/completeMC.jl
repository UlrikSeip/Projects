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
    return A*sin(omega*t) + a0
end

function getf(t, f0, F, omegaf)
    F*sin(omegaf*t) + f0
end

function MC_sol(a0, A, omega, b, c, d, di, bi, f0, ft0, ft1, df, F, omegaf, S0, I0, R0, T, sims, filename, dt = 0, antivac = false, camp = false)
    """
    Takes the inputs rate of transmission a0, variance in transmission A, frequency of transmission function omega,
    rate of recovery b, rate of immunity loss c, rate of death d, rate of death by sickness di, rate of birth bi,
    rate of vaccination f0, variance in vaccination F, frequency of variance in vaccination omegaF,
    initial values of susceptible , infected and recovered individuals, S0, I0 and R0 respectively,
    the total simulation time T in days, and file filename. Can allso take a custom dt.
    Plots the resulting arrays and saves the file as filename.
    """
    a = geta(0, a0, A, omega)
    N0 = S0+I0+R0 #the total population
    tmins = [4/(a*N0), 1/(b*N0), 1/c*N0]
    dt == 0 ? dt = minimum(tmins) : dt = dt
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
            f = getf(t[q], f0+(t[q]*df), F, omegaf)
            if camp
                f = 0
                if t[q] > ft0 && t[q] < ft1
                    f = f0
                end
            end
            if antivac
                if I[q] < 5
                    f = 0
                    f0 = 0
                    df = 0
                end
            end
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
    is = (-(b/a)+1)/(1 + (b/c))
    rs = (b/c)*(1-(b/a))/(1 + (b/c))


    println("Final tally: N=$(sum(Savg[end]+Iavg[end]+Ravg[end]))")
    println("Expected: S=$(b/a), I=$is, R=$rs")
    println("Found: S=$(Savg[end]/Navg[end]), I=$(Iavg[end]/Navg[end]), R=$(Ravg[end]/Navg[end])")
    println("Number: S=$(Savg[end]), I=$(Iavg[end]), R=$(Ravg[end])")
    println("Standard deviations: S=$Sstd, I=$Istd, R=$Rstd")
    println()
    plt.figure()
    plt.plot(t,Savg)
    plt.plot(t,Iavg)
    plt.plot(t,Ravg)
    plt.grid()
    plt.legend(["Susceptible", "Infected", "Recovered", "Total"]) #"Dead",
    plt.xlabel("Days")
    plt.ylabel("People")
    plt.savefig(filename)
    #plt.show()

end

a0 = 4                         #lingerling chem trails
A = 2                          #A = 0 to disable oscillator (size of chem trails)
omega = 2*pi/365.25            #frequency of chem-trails
b = 1                          #effectiveness of healing crystals
c = 0.5                        #levels of atheism
d = 0.00002242299              #government abductions
di = 0                       #rate of purge
bi = 0.00002948891             #newborn heathens
f0 = 0                       #f = 0 for effective anti-vac campaigns
df = 0                         
ft0 = 2
ft1 = 9
F = 0                          
omegaf = 2*pi/365.25
S0 = 300            
I0 = 100
R0 = 0
T = 356.25*2
sims = 1000
dt = 0
antivac = false
camp = false
filename = "opp_d_C_MC.pdf"

MC_sol(a0, A, omega, b, c, d, di, bi, f0, ft0, ft1, df, F, omegaf, S0, I0, R0, T, sims, filename, dt, antivac, camp)
