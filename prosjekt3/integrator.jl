using PyCall
@pyimport matplotlib.pyplot as plt
import LinearAlgebra: norm

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


function forward_euler(vel0, pos0, t, dt, endvalue)
    #NB!! simplified for circular motion in G-field
    #vel0 and pos0 should be (x, y, z) arrays with initial values for vel an pos
    #endvalues is a bool, return only resulting value after time t
    #otherwise return entire pos and vel array
    length = int(t/dt)
    vel = ones(3, length)
    vel[:, 1] = vel0
    pos = ones(3, length)
    pos[:, 1] = pos0
    for i = 2:length
        rip = norm(pos[:, i-1]) #r[i-1]
        th = arccos(pos[:, i-1]/rip)
        v[:, i] = vel[:, i-1]-4*pi*dt/(rip^2))*cos(th)
        pos[:, i] = pos[:, i-1] + dt*v[:, i]
    end
    if endvalue
        return pos[:, -1], vel[:, -1]
    end
    return pos, vel
end

function velocity_verlet(vel0, pos0, dvdt, t, dt, endvalue)
    #vel0 and pos0 should be (x, y, z) arrays with initial values for vel an pos
    #endvalues is a bool, return only resulting value after time t
    #otherwise return entire pos and vel array
    length = int(t/dt)
    vel = ones(3, length)
    vel[:, 1] = vel0
    pos = ones(3, length)
    pos[:, 1] = pos0
    for i = 2:length
        rip = norm(pos[:, i-1]) #r[i-1]
        th = arccos(pos[:, i-1]/rip)
        aip = -4*pi^2/(rip^2)*cos(th)   #a[i-1]
        pos[:, i] = pos[:, i-1] + dt*vel[:, i-1] + ((dt^2)/2)*aip
        ri = rip = norm(pos[:, i])
        ai = -4*pi^2/(ri^2)*cos(th)
        vel[:, i] = vel[:, i-1] + dt*(ai + aip)/2
    end
    if endvalue
        return pos[:, -1], vel[:, -1]
    end
    return pos, vel
end

#read info from file
plt.#3dplot somehow
plt.show()