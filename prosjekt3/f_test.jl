include("integrator.jl")
using PyCall
using PyPlot
using LinearAlgebra
#using Math
@pyimport matplotlib.pyplot as plotter
@pyimport numpy as np
import LinearAlgebra: norm

function acc_fs(pos, par)
    """
    Finds the acceleration for mulitple planets when the sun is fixed in the centre
    """
    #some constants
    s_mass = 1
    masses = par
    G = -4*(pi^2)
    G0 = -6.67e-11
    nr_pl = Int64(length(pos)/3) #the number of planets
    pl = zeros(Int64(nr_pl), 3) #an arrray to better hold all the positions of the planets

    #inserts the positions of the planets in pl
    i = 3
    while i < length(pos)+1
        temp = pos[i-2:i]
        pl[Int64(i/3),:] = temp
        i+=3
    end
    a = zeros(length(pos)) #an array to hold all the accelerations in the same style as pos
    ting = true
    #goes through all the planets
    for i = 1:nr_pl
        k = Int64(3*i)
        #finds the acceleration from the sun
        r = norm(pl[i,:])
        th = pl[i,:]/r
        a[k-2:k] = G*s_mass*th/((r)^2)

        #finds the acceleration from the other planets
        for j = 1:nr_pl
            if i != j
                dis = pl[i,:] - pl[j,:]
                r = norm(dis)
                th = dis/r
                a[k-2:k] += G*parse(Float64,masses[j])*th/((r)^2)
            end
        end
    end
    return a
end

function acc_nfs(pos, par)
    """
    Finds the acceleration for mulitple planets when the sun is not fixed in the centre, but treated as an object
    """
    #some constants
    masses = par
    G = -4*(pi^2)
    nr_pl = Int64(length(pos)/3) #the number of objects
    pl = zeros(Int64(nr_pl), 3) #an arrray to better hold all the positions of the objects

    #inserts the positions of the objects in pl
    i = 3
    while i < length(pos)+1
        temp = pos[i-2:i]
        pl[Int64(i/3),:] = temp
        i+=3
    end

    a = zeros(length(pos)) #an array to hold all the accelerations in the same style as pos
    #goes through all the objects
    for i = 1:nr_pl
        k = Int64(3*i)  #dimensions? ------------------------

        #finds the acceleration from the other objects
        for j = 1:nr_pl
            if i != j
                dis = pl[i,:] - pl[j,:]
                r = norm(dis) #-sjekk om denne er riktig -----------------------------
                th = dis/r
                a[k-2:k] += G*parse(Float64,masses[j])*th/((r)^2) #y not a[1:k]? --------------------------
                #println(float(masses[j]))
            end
            #print a and check if reasonable ------------------------
        end
    end
    return a
end


"""

earthM = 5.97219e24
earthPos0 = [9.413801075750535E-01, 3.379019986046322E-01, -9.334104672733438E-05]
earthVel0 = [-5.994522787486753E-03, 1.617377250092178E-02, -1.732657683299539E-07]

saturnM = 5.6834e26
saturnPos0 = [1.543954489897663E+00, -9.936037960121499E+00, 1.113028547106575E-01]
saturnVel0 = [5.205269695259624E-03, 8.394344861930620E-04, -2.220728858861169E-04]

sunM = 2e30

pos0 = zeros(6)
pos0[1:3] = earthPos0
pos0[4:6] = saturnPos0
println(pos0)

vel0 = zeros(6)
vel0[1:3] = earthVel0
vel0[4:6] = saturnVel0


a = acc_fs(pos0, [earthM, saturnM])
println(a)


stopTime = 2
pos, vel = velocity_verlet(pos0, vel0, stopTime, stopTime/1e5, acc_fs, [sunM, [earthM, saturnM]])
println(size(pos))


"""