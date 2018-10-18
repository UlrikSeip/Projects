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


function forward_euler(vel0, pos0, t, dt, func, endvalue = false) 
"""
NB!! simplified for circular motion in G-field
vel0 and pos0 should be (x, y, z) arrays with initial values for vel an pos
endvalues is a bool, return only resulting value after time t
otherwise return entire pos and vel array
"""
    
    #creating initial values
    length = Int64(round(t/dt))
    vel = ones(3, length)
    vel[:, 1] = vel0
    pos = ones(3, length)
    pos[:, 1] = pos0

    #integration loop
    for i = 2:length
        a = func(pos[:, i-1])
        vel[:, i] = vel[:, i-1] + a*dt #new vel
        pos[:, i] = pos[:, i-1] + dt*vel[:, i] #new pos
    end
    #return related stuff
    if endvalue
        return pos[:, end], vel[:, end]
    end
    return pos, vel
end

function velocity_verlet(vel0, pos0, t, dt, func, endvalue = false)
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
    ai = func(pos[:, 1]) 

    #integration loop
    for i = 2:length
        aip = ai   #current acceleration
        pos[:, i] = pos[:, i-1] + dt*vel[:, i-1] + ((dt^2)/2)*aip  #new position, don't know why, but we need to write [aip[1], aip[2], aip[3]] instead of aip
        ai = func(pos[:, i]) #new acceleration based on new radius
        vel[:, i] = vel[:, i-1] + dt*(ai + aip)/2   #new velocity based on new acceleration

    end
    #return related stuff
    if endvalue
        return pos[:, end], vel[:, end]
    end
    return pos, vel
end

function filewriter(dataArray, filename = "../prosjekt3/orbits.txt")
    #dataArray should be a 3d array, and is written to filename as a[:, 1]\n, a[:, 2]\n...
    f = open(filename,"w")
    for i = 1:Int64(length(dataArray)/3)
        writedlm(f, dataArray[:, i])
    end
    close(f)
end

