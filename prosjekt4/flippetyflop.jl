println("Importing PyCall")
@time using PyCall
println("Importing PyPlot")
@time using PyPlot
const plt = PyPlot
println("Importing small stuff")
@time using Random
@time using Printf
@time using Dates
@time import LinearAlgebra: norm
@time import Statistics
@time import Base.Threads
@pyimport matplotlib.animation as anim
const kb = 8.6173303E-5

println()

const threads = Sys.CPU_THREADS
Threads.nthreads() = threads

function laticemaker(dim1, dim2, seed = 0, justones = false)
    #Generates a random lattice with option of a seed. Seed = 0 --> no seed
    seed != 0 ? Random.seed!(seed) : 42
    justones ? latice = latice = ones((dim1, dim2)) : rand([-1, 1], dim1, dim2)
end

function get_E(M, J = 1) 
    #=
    Sums total energy over all particles
    M = lattice
    J = constant
    =#

    L = length(M[1, :])
    E_tot = 0
    @inbounds for i = 1:2:L
        @inbounds for j = 1:2:L
            ns = zeros(4)   #left, top, right, bottom
            #bool? true:false
            i > 1 ? ns[1] = M[i-1, j] : ns[1] = M[L, j]
            j > 1 ? ns[2] = M[i, j-1] : ns[2] = M[i, L]
            i < L ? ns[3] = M[i+1, j] : ns[3] = M[1, j]
            j < L ? ns[4] = M[i, j+1] : ns[4] = M[i, 1]

            E_tot += sum(ns*M[i, j])
        end
    end
    @inbounds for i = 2:2:L
        @inbounds for j = 2:2:L
            ns = zeros(4)   #left, top, right, bottom
            #bool? true:false
            i > 1 ? ns[1] = M[i-1, j] : ns[1] = M[L, j]
            j > 1 ? ns[2] = M[i, j-1] : ns[2] = M[i, L]
            i < L ? ns[3] = M[i+1, j] : ns[3] = M[1, j]
            j < L ? ns[4] = M[i, j+1] : ns[4] = M[i, 1]

            E_tot += sum(ns*M[i, j])
        end
    end
    return -J*Int64(E_tot)
end

function get_dE(s, ns, J = 1)
    #calculates around a single spin s and its neigbours ns [left, top, right, bottom]
    return(Int64(2J*sum(s*ns)))
end

function get_M(M)
    #Returns total magnetic moment of matrix M
    M_tot = 0
    L = length(M[:, 1])
    for i = 1:L
        for j = 1:L
            M_tot += M[i,j]
        end
    end
    return M_tot
end

function flip(M, flips)
    #=
    Changes sign of requested elements in M
    M     = matrix to flip
    flips = array [[xs], [ys]] to flip
    =#
    for i = 1:length(flips[:, 1])
        x = flips[i, 1]
        y = flips[i, 2]
        M[x, y] *= -1
    end
    return M
end

function rand_ind(N, L)
    #generates N random pairs of indices in 1:L
    indices = zeros(Int64, (N, 2))
    for i = 1:N
        indices[i, :] = rand(1:L, 2)
    end
    return indices
end

function print_an_2x2(T, J = 1)
    #prints analytical values for the 2x2 lattice
    Z = 4*(cosh(8J/T)+3)
    E = (32/Z)*sinh(-8J/T)
    Ma = (8/Z)*(exp(8J/T) + 2)
    Cv = (kb/(T^2))*((128/Z)*cosh(-8J/T) - E^2)
    chi = (32/(Z*T))*(exp(8J/T) + 1 - (2/Z)*(exp(16J/T) + 2*exp(8J/T) + 4))
    println("2x2 Analytical Data")
    println("⟨E⟩, ⟨|M|⟩, Cᵥ, χ")
    println("$E, $Ma, $Cv, $chi")
end

function test_get_E()
    #Test the "get_E" function
    M = [-1 1; 1 -1]
    expected = 8
    calculated = get_E(M)
    assert(expected == calculated, "get_E: expected $expected != computed $calculated")
end

function test_get_M()
    #Test the "get_M" function
    M = [1 1 1 1; -1 -1 -1 -1; 1 -1 1 -1; -1 1 1 1]
    expected = 2
    calculated = get_M(M)
    assert(expected == calculated, "get_M: expected $expected != computed $calculated")
end

function test_flip()
    #Test the "flip" function
    M = [1 1 1 1; -1 -1 -1 -1; 1 -1 1 -1; -1 1 1 -1]
    toflip = [1 1; 1 2; 1 4; 3 2]
    expected = [-1 -1 1 -1; -1 -1 -1 -1; 1 1 1 -1; -1 1 1 -1]
    calculated = flip(M, toflip)
    assert(expected == calculated, "flip: expected $expected != computed $calculated")
end

function test_all()
    #Runs all test functions
    tests = [test_get_E, test_flip, test_get_M]
    N = length(tests)
    for i = 1:N
        println("Testing [$i/$N]")
        tests[i]()
    end
    println()
end

function assert(toTest, name)
    #Works sortof like the assert method in Python
    if toTest == false
        println()
        error("Error in function: \"$name\"")
    end
end

function print_matrix(M)
    #Neatly prints matrix M
    N = size(M)[1]
    println("------Matrix------")
    for i=1:N
        for j=1:N
            print("$(@sprintf("%7.2f ", M[i,j]))")
        end
        println()
    end
    println()
end

function drawLattice(M)
    #Shows lattice M
    plt.imshow(M, cmap = "binary")
    plt.show()
end

function plot_sim(M, L, t0, t1, N, J = 1)
    #shows an image of the final lattice
    @time energy, magmom = sim(M, L, N, t0, t1, J)
    cycles = collect(0:N)
    plt.plot(cycles, energy)
    plt.xlabel("Cycles")
    plt.ylabel("Energy")
    plt.show()
    plt.plot(cycles, magmom)
    plt.xlabel("Cycles")
    plt.ylabel("Magnetic moment")
    plt.show()
end

function sim(M, L, N, T, J = 1)
    #simple simulation, doesnt return anything but the basic values
    Es = zeros(N+1)
    Es[1] = get_E(M, J)
    Ms = zeros(N+1)
    Ms[1] = get_M(M)
    for i = 1:N
        ind = rand_ind(L^2, L)
        M = step(M, T, ind, J)
        Es[i+1] = get_E(M, J)
        Ms[i+1] = get_M(M)
    end
    return Es, Ms, M
end


function step(M, t, ind, J = 1)
    #one step of the algorithm
    L = length(M[1, :])
    final_ind = zeros(Int64, size(ind))
    ficount = 0
    w = zeros(17)
    for i=-8:4:8
        w[i+9] = exp(-J*i/t)
    end
    @inbounds for k = 1:size(ind)[1]
        i, j = ind[k, :]
        ns = zeros(4)   #left, top, right, bottom
        i > 1 ? ns[1] = M[i-1, j] : ns[1] = M[L, j]
        j > 1 ? ns[2] = M[i, j-1] : ns[2] = M[i, L]
        i < L ? ns[3] = M[i+1, j] : ns[3] = M[1, j]
        j < L ? ns[4] = M[i, j+1] : ns[4] = M[i, 1]

        dE = get_dE(M[i, j], ns, J)
        if (dE < 0) || (rand() < w[Int64(dE+9)])
            ficount += 1
            final_ind[ficount, :] = [i, j]
        end
    end
    return flip(M, final_ind[1:ficount, :])
end

function get_data(L, T, steps, J = 1)
    #one run
    M = laticemaker(L, L)
    E = get_E(M)
    Ma = get_M(M)
    Etot = 0
    Matot = 0
    EtotSq = 0
    MatotSq = 0
    w = zeros(17)
    for i=-8:4:8
        w[i+9] = exp(-J*i/T)
    end
    @inbounds for q = 1:steps
        ind = rand_ind(L^2, L)
        Threads.@threads for k = 1:size(ind)[1]
            i, j = ind[k, :]
            ns = zeros(4)   #left, top, right, bottom
            i > 1 ? ns[1] = M[i-1, j] : ns[1] = M[L, j]
            j > 1 ? ns[2] = M[i, j-1] : ns[2] = M[i, L]
            i < L ? ns[3] = M[i+1, j] : ns[3] = M[1, j]
            j < L ? ns[4] = M[i, j+1] : ns[4] = M[i, 1]

            dE = get_dE(M[i, j], ns, J)
            if (dE < 0) || (rand() < w[Int64(dE+9)])
                E += dE
                M[i, j] *= -1
                Ma += 2*M[i, j]
            end
        end
        Etot += E
        EtotSq += E^2
        Matot += abs(Ma)
        MatotSq += Ma^2
    end
    Cv = EtotSq/steps - Etot^2/(steps^2)
    Chi = MatotSq/steps - Matot^2/(steps^2)
    return Etot/steps, Matot/steps, Cv, Chi
end


function get_values(L, T0, T1, dt, steps, runs, J = 1)
    #=
    Simulates, and finds average values for "runs" runs
    =#   
    Ts = T0:dt:T1
    println("Initiating main simulation at: ", Dates.format(now(), "HH:MM:SS"))
    println("Run info: $L x $L, T($T0-$T1), dt = $(round(dt; digits = 5))")
    p = 10
    fudgefactor = 1
    t = @elapsed Threads.@threads for i = 1:p
        get_data(L, Ts[1], steps, J)
    end
    totsteps = length(Ts)*runs
    println("Estimated seconds until completion: $((t/p)*(totsteps*fudgefactor))")
    Es = zeros(length(Ts))
    Mas = zeros(length(Ts))
    Cvs = zeros(length(Ts))
    Chis = zeros(length(Ts))
    @inbounds for i = 1:length(Ts)
        vals = zeros(4)
        Threads.@threads for q = 1:runs
            vals = vals.+get_data(L, Ts[i], steps, J)
        end
        vals = vals./runs
        Es[i] = vals[1]
        Mas[i] = vals[2]
        Cvs[i] = vals[3]
        Chis[i] = vals[4]
    end
    println()
    return Es, Mas, Cvs, Chis, Ts
end

function gifSim(savefile, M, L, T0, T1, N, J = 1)
    #create a mp4 file of a simulation.
    skip = 1
    println("L $L steps $N t0 $T0 t1 $T1.mp4")
    fig = figure()
    im =  imshow(M, cmap = "binary")
    ims = [PyCall.PyObject[im]]
    dT = (T1-T0)/N
    for i in 1:N-1
        ind = rand_ind(L^2, L)
        M = step(M, T0 + i*dT, ind, J)
        im = imshow(M, cmap = "binary")
        push!(ims, PyCall.PyObject[im])
    end
    ani = anim.ArtistAnimation(fig, ims, interval=50, blit=true)
    ani[:save](savefile, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"]);
    println()
end

function gifSims(Ls, T0s, T1s, stepss)
    #run gifSim with several different combinaitons of variables
    for steps in stepss
        for T0 in T0s
            for T1 in T1s
                for L in Ls
                    M = laticemaker(L, L)
                    savefile = "L $L steps $steps T0 $T0 T1 $T1.mp4"
                    @time gifSim(savefile, M, L, T0, T1, steps)
                end
            end
        end
    end
end

function plot_values(Es, Mas, Cvs, Chis, Ts, Ls, steps, runs, J = 1, analytical = false)
    #=
    Plots 2d arrays of the first 4 arguments against T. Saves as png.
    If analytical = true, overlays analytical solution for 2x2
    =#
    T0 = minimum(Ts)
    T1 = maximum(Ts)

    Esa = zeros(length(Ts))
    Msa = zeros(length(Ts))
    Cvsa = zeros(length(Ts))
    Chisa = zeros(length(Ts))

    if analytical
        for i = 1:length(Ts)
            T = Ts[i]
            Z = 4*(cosh(8J/T)+3)
            Esa[i] = (128/Z)*(exp(8J/T) + exp(-8J/T))
            Msa[i] = (32/Z)*(exp(8J/T) + 1)
            Esaavg = -(16/Z)*(exp(8J/T) - exp(-8J/T))
            Masaavg = (8/Z)*(exp(8J/T) + 2)
            Cvsa[i] = (Esa[i] - Esaavg^2)/(T^2)
            Chisa[i] = (Msa[i] - Masaavg^2)/(T^2)
        end
    end
    println(size(Ts))
    println(size(Es))
    for i = 1:length(Es[:, 1])
        plt.plot(Ts, Es[i, :], label = "L = $(Ls[i])")
    end
    plt.xlabel("Temperature [kB*K]")
    plt.ylabel("Mean energy")
    plt.yscale("log")
    plt.title("\$ \\langle E \\rangle \$")
    plt.legend()
    plt.savefig("E_$(Ls[1])_to_$(Ls[end])_$(steps)_T($(T0)_$(T1)).png")
    plt.close()

    for i = 1:length(Mas[:, 1])
        plt.plot(Ts, Mas[i, :], label = "L = $(Ls[i])")
    end
    plt.xlabel("Temperature [kB*K]")
    plt.ylabel("Mean energy")
    plt.yscale("log")
    plt.title("\$ \\langle |M| \\rangle \$")
    plt.legend()
    plt.savefig("Ma_$(Ls[1])_to_$(Ls[end])_$(steps)_T($(T0)_$(T1)).png")
    plt.close()

    for i = 1:length(Cvs[:, 1])
        plt.plot(Ts, Cvs[i, :], label = "L = $(Ls[i])")
    end
    plt.xlabel("Temperature [kB*K]")
    plt.ylabel("Mean energy")
    plt.yscale("log")
    plt.title("\$ C_V \$")
    plt.legend()
    plt.savefig("Cv_$(Ls[1])_to_$(Ls[end])_$(steps)_T($(T0)_$(T1)).png")
    plt.close()
    
    for i = 1:length(Chis[:, 1])
        plt.plot(Ts, Chis[i, :], label = "L = $(Ls[i])")
    end
    plt.xlabel("Temperature [kB*K]")
    plt.ylabel("Mean energy")
    plt.yscale("log")
    plt.title("\$ \\chi \$")
    plt.legend()
    plt.savefig("Chi_$(Ls[1])_to_$(Ls[end])_$(steps)_T($(T0)_$(T1)).png")
    plt.close()
end

function print_values(Es, Mas, T, steps)
    #prints analytical values 
    Z = 4*(cosh(8J/T)+3)
    Esa= (128/Z)*(exp(8J/T) + exp(-8J/T))
    Msa= (32/Z)*(exp(8J/T) + 1)
    Esaavg = -(16/Z)*(exp(8J/T) - exp(-8J/T))
    Masaavg = (8/Z)*(exp(8J/T) + 2)
    Cvsa= (Esa - Esaavg^2)/(T^2)
    Chisa = (Msa - Masaavg^2)/(T^2)
    
end


function main()
    test_all()
    Ls = [20, 40, 60, 80, 100]
    steps = 800
    runs = 5
    J = 1
    lenTs = 30
    Es = zeros((length(Ls), lenTs+1))
    Mas = zeros((length(Ls), lenTs+1))
    Cvs = zeros((length(Ls), lenTs+1))
    Chis = zeros((length(Ls), lenTs+1))
    Ts = 0
    T0 = 1
    print_an_2x2(T0)
    println(get_values(2, 1, 2, 1, 1000, 10000, J))
    println()
    for i = 1:length(Ls)
        L = Ls[i]
        T0 = 1
        T1 = 3
        dt = (T1-T0)/lenTs
        Es[i, :], Mas[i, :], Cvs[i, :], Chis[i, :], Ts = get_values(L, T0, T1, dt, steps, runs, J)
    end
    plot_values(Es, Mas, Cvs, Chis, Ts, Ls, steps, runs, J)

    """
    @inbounds for L in Ls
        T0 = 1
        T1 = 3
        dt = (T1-T0)/100
        Es, Ms, Cvs, Chis, Ts = get_values(L, T0, T1, dt, steps, runs, J)
        plot_values(Es, Ms, Cvs, Chis, Ts, L, steps, runs, J)
        
        T0 = 2
        T1 = 2.3
        dt = (T1-T0)/100
        Es, Ms, Cvs, Chis, Ts = get_values(L, T0, T1, dt, steps, runs, J)
        plot_values(Es, Ms, Cvs, Chis, Ts, L, steps, runs, J, anal)
    end
    """
    #T0 = 1
    #T1 = 1
    #gifSims(Ls, T0s, T1s, stepss)
    #Es, Ms, Cvs, Chis, Ts = get_values(L, T0, T1, dt, steps, runs, J)
end

main()
