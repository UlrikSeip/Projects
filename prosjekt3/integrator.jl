using PyCall
using PyPlot
using LinearAlgebra
using DelimitedFiles
#using Math
@pyimport matplotlib.pyplot as plotter
@pyimport numpy as np
import LinearAlgebra: norm

#by starting a julia shell, and typing "]" you cna then "add PyCall"
#and "add PyPlot" to install packages

#julia integrator class, probably, eventually
#essentially an ode solver with several different integration methods
#object oriented syntax in julia:
"""
type MyType
a::Int64
b::Float64
end

x = MyType(3, 4)

x.a

---------------------
function double(x::MyType)
    x.a *= 2
end
"""

function acc_sirc(pos)
    """
    Finds the acceleration for an object in a constant circular motion, in a spesific position.
    Takes the position of the object.
    """    
    r = norm(pos)
    th = pos/r
    a = 4*pi^2/(r^2)*th
    return a
end


function forward_euler(vel0, pos0, t, dt, endvalue = false)
"""
NB!! simplified for circular motion in G-field
vel0 and pos0 should be (x, y, z) arrays with initial values for vel an pos
endvalues is a bool, return only resulting value after time t
otherwise return entire pos and vel array
"""
    
    #creating initial values
    length = Int64(t/dt)
    vel = ones(3, length)
    vel[:, 1] = vel0
    pos = ones(3, length)
    pos[:, 1] = pos0

    #integration loop
    for i = 2:length
        #rip = norm(pos[:, i-1]) #new r
        #print(pos[:, i-1])
        #th = pos[:, i-1]/rip #new angle
        a = acc_sirc(pos[:, i-1])
        vel[:, i] = vel[:, i-1]-a*dt #new vel
        pos[:, i] = pos[:, i-1] + dt*vel[:, i] #new pos
    end
    #return related stuff
    if endvalue
        return pos[:, end], vel[:, end]
    end
    return pos, vel
end

function velocity_verlet(vel0, pos0, t, dt, endvalue = false)
"""
vel0 and pos0 should be (x, y, z) arrays with initial values for vel an pos
endvalues is a bool, return only resulting value after time t
otherwise return entire pos and vel array
"""

    #initial values
    length = Int64(round(t/dt))
    vel = ones(3, length)
    vel[:, 1] = vel0
    pos = ones(3, length)
    pos[:, 1] = pos0
    #rip = norm(pos[:, 1])
    ai = acc_sirc(pos[:, 1]) 

    #integration loop
    for i = 2:length
        #rip = norm(pos[:, i-1]) #current radius
        #th = pos[:, i-1]/rip #current angle

        aip = ai   #current acceleration
        #println(pos[:, i-1], vel[:, i-1], aip)
        pos[:, i] = pos[:, i-1] + dt*vel[:, i-1] + ((dt^2)/2)*[aip[1], aip[2], aip[3]]  #new position, don't know why, but we need to write [aip[1], aip[2], aip[3]] instead of aip
        #if i == 2
        #    println(th)
        #    println(aip)
        #    println(pos[:, i])
        #end
        #rip = norm(pos[:, i])  #new radius
        ai = acc_sirc(pos[:, i]) #new acceleration based on new radius
        vel[:, i] = vel[:, i-1] + dt*([ai[1], ai[2], ai[3]] + [aip[1], aip[2], aip[3]])/2   #new velocity based on new acceleration
    end
    #return related stuff
    if endvalue
        return pos[:, end], vel[:, end]
    end
    return pos, vel
end

function filewriter(dataArray, filename = "orbits.txt")
    #dataArray should be a 3d array, and is written to filename as a[:, 1]\n, a[:, 2]\n...
    f = open(filename,"w")
    for i = 1:Int64(length(dataArray)/3)
        writedlm(f, dataArray[:, i])
    end
    close(f)
end

#i think we need a function for acceleration. our orbits mainly have positive values...

"""
function aFunk(r, ang, mass = 1)
G = 0.01720209895
return G*mass/(r^2)*cos.(ang)
end
"""



v0 = [2.24e-03, 1.51e-02, 2.61e-0]
x0 = [9.41e-01, 3.38e-01, -9.33e-05]
stopTime = 3
#3.154e+7
pos, vel = forward_euler(v0, x0, stopTime, stopTime/36000)
filewriter(pos)
#read info from file

println(size(pos))
#println(pos)
#plotter.plot(pos[1], pos[2]) #np.linspace(0, stopTime, length(pos[1])))
#when written to file, the data is better presented by the plitting function in solarsystem.py
#plotter.show()
