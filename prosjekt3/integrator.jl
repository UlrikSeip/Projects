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

#planetinfo
sunM = 1
#sunM = 2e30
mercuryM = 3.302e23
mercuryPos0 = [-1.741425018981730E-01, -4.236575369435728E-01, -1.932082627282657E-02]
mercyryVel0 = [2.042279845238355E-02, -9.147939192029718E-03, -2.621780116213209E-03]
venusM = 48.685e23
venusPos0 = [7.108874104850517E-01, 1.475289827913959E-01, -3.917897844928585E-02]
venusVel0 =[-3.997930253479682E-03, 1.975505990962758E-02, 5.014443552290667E-04]
earthM = 5.97219e24
earthPos0 = [9.413801075750535E-01, 3.379019986046322E-01, -9.334104672733438E-05]
earthVel0 = [-5.994522787486753E-03, 1.617377250092178E-02, -1.732657683299539E-07]
marsM = 6.4171e23
marsPos0 = [1.375357774690336E+00, -1.627517936385902E-01, -3.738675132962930E-02]
marsVel0 = [2.243217631343320E-03, 1.508628660091320E-02, 2.610262676274213E-04]
jupiterM = 1898.13e24
jupiterPos0 = [-2.666952709077877E+00, -4.655671225645230E+00, 7.896515774211305E-02]
jupiterVel0 = [6.458958874387921E-03, -3.390642961368397E-03, -1.303431975919576E-04]
saturnM = 5.6834e26
saturnPos0 = [1.543954489897663E+00, -9.936037960121499E+00, 1.113028547106575E-01]
saturnVel0 = [5.205269695259624E-03, 8.394344861930620E-04, -2.220728858861169E-04]
uranusM = 86.813e24
uranusPos0 = [1.717792200893557E+01, 9.994448097964824E+00, -1.854226138335391E-01]
uranusVel0 = [-2.006740577381741E-03, 3.216232224304096E-03, 3.798323252188187E-05]
neptuneM = 102.413e24
neptunePos0 = [2.891951044185697E+01, -7.725591586322342E+00, -5.073855242738878E-01]
neptuneVel0 = [7.889374352298129E-04, 3.051510114439145E-03, -8.066487877030014E-05]
plutoM = 1.307e22
plutoPos0 = [1.164087973998699E+01, -3.157554534292389E+01, 1.155639944820950E-02]
plutoVel0 = [3.024320808138397E-03, 4.313236772252202E-04, -9.237227795187408E-04]



function acc_sirc(pos)
    """
    Finds the acceleration for an object in a constant circular motion, in a spesific position.
    Takes the position of the object.
    """    
    r = norm(pos)
    th = pos/r
    a = -4*pi^2/(r^2)*th
    return a
end


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

function filewriter(dataArray, filename = "orbits.txt")
    #dataArray should be a 3d array, and is written to filename as a[:, 1]\n, a[:, 2]\n...
    f = open(filename,"w")
    for i = 1:Int64(length(dataArray)/3)
        writedlm(f, dataArray[:, i])
    end
    close(f)
end

#i think we need a function for acceleration. our orbits mainly have positive values...


function aFunk(pos, mass = sunM)
    #G = 6.67408e-11
    G = -4*(pi^2)
    #auToM = 149.60*10^9
    r = norm(pos)
    th = pos/r
    a = G*mass*th/((r)^2)
    return a
end




x0 = [9.41e-01, 3.38e-01, -9.33e-05]    #earth
v0 = 365.2422*[-5.99e-03, 1.62e-02, -1.73e-07]    #earth
#x0 = [1,0,0]
#v0 = [0,2,0]
stopTime = 10    #nr of years or days for the simulation 
#3.154e+7
pos, vel = velocity_verlet(365.2422*earthVel0, earthPos0, stopTime, stopTime/1e5, aFunk)
filewriter(pos)
#read info from file

println(size(pos))
#println(pos)
#plotter.plot(pos[1], pos[2]) #np.linspace(0, stopTime, length(pos[1])))
#when written to file, the data is better presented by the plitting function in solarsystem.py
#plotter.show()