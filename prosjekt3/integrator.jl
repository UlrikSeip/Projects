using DelimitedFiles
include("afunks.jl")
#by starting a julia shell, and typing "]" you cna then "add PyCall"
#and "add PyPlot" to install packages

function forward_euler(vel0, pos0, t, dt, func, par, endvalue = false) 
"""
vel0 and pos0 should be (x, y, z) arrays with initial values for vel an pos
endvalues is a bool, return only resulting value after time t
otherwise return entire pos and vel array
"""
    
    #creating initial values
    len = Int64(round(t/dt))
    vel = ones(3, len)
    vel[:, 1] = vel0
    pos = ones(3, len)
    pos[:, 1] = pos0

    #integration loop
    for i = 2:len
        a = func(pos[:, i-1], par)
        vel[:, i] = vel[:, i-1] + a*dt #new vel
        pos[:, i] = pos[:, i-1] + dt*vel[:, i] #new pos
    end
    #return related stuff
    if endvalue
        return pos[:, end], vel[:, end]
    end
    return pos, vel
end

function velocity_verlet(vel0, pos0, t, dt, func, par, endvalue = false)
"""
vel0 and pos0 should be [[planets],[x, y, z]] arrays with initial values for vel an pos
endvalues is a bool, return only resulting value after time t
otherwise return entire pos and vel array
"""

    #initial values
    len = Int64(round(t/dt))
    pl_len = length(vel0)
    vel = ones(pl_len, len)
    vel[:, 1] = vel0
    pos = ones(pl_len, len)
    pos[:, 1] = pos0
    println(pos0)
    #rip = norm(pos[:, 1])
    #println(size(pos))
    ai = func(pos[:, 1], par) 
    println(ai)

    #integration loop
    for i = 2:len
        aip = ai   #current acceleration
        #println(size(aip))
        #println(size(pos))
        pos[:, i] = pos[:, i-1] + dt*vel[:, i-1] + ((dt^2)/2)*aip  #new position
        ai = func(pos[:, i], par) #new acceleration based on new radius
        vel[:, i] = vel[:, i-1] + dt*(ai + aip)/2   #new velocity based on new acceleration

    end
    #return related stuff
    if endvalue
        return pos[:, end], vel[:, end]
    end
    return pos, vel
end

#6191060/160

function filewriter(dataArray, filename = "../prosjekt3/orbits.txt")
    #dataArray should be a 3d array, and is written to filename as a[:, 1]\n, a[:, 2]\n...
    f = open(filename,"w")
    for i = 1:Int64(length(dataArray)/3)
        writedlm(f, dataArray[:, i])
    end
    close(f)
end
