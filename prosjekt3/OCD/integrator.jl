include("afunks.jl")
include("plotter.jl")
import PyPlot
const plt = PyPlot
using ArgParse
using DelimitedFiles
using NPZ
#using PyPlot


"""
Simply run the file with the required command line arguments to integrate the supplied values.
The values should be stored in a file with the format:
"""


function parse_commandline() #equivalent to argparse in python
    args = ArgParseSettings()

    @add_arg_table args begin
        "readfile"
            help = "name of file with initial info with alternating x0 and v0 for several objects"
            required = true
            arg_type = String
        "writefile"
            help = "name of destination file for integrated values"
            required = true
            arg_type = String
        "t"
            help = "time for simulation"
            required = true
            arg_type = Float64
        "dt"
            help = "dt"
            required = true
            arg_type = Float64
        "plot"
            help = "plot function"
            required = false
            default = false
            arg_type = Bool
        "masses"
            help = "array with masses of planets in the order they appear in readfile"
            required = false
            default = [3.302e23, 48.685e23, 5.97219e24, 6.4171e23, 1.89813e27, 5.6834e26, 86.813e24, 102.413e24, 1.307e22]
            #arg_type = Array
            #det mangler én masse. Dunno hvilket
        "names"
            help = "array with names of planets in the order they appear in readfile"
            required = false
            default = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptun", "Pluto"]
            #arg_type = Array
        "acc_func"
            help = "a number representing which acceleration function you want to use; 1=aFunk, 2=acc_fs, 3=acc_nfs, 4=moreBodyFunc"
            required = false
            default = 1
            arg_type = Int64
        "int_func"
            help = "integration method"
            required = false
            default = 1
            arg_type = Int64
        "plot_dim"
            help = "plot in 2 or 3 dimensions"
            required = false
            default = 2
    end
    return parse_args(args)
end

function parse() #reads and returns initial info from file 
    parsed_args = parse_commandline()
    data = npzread(parsed_args["readfile"])
    writefile = parsed_args["writefile"]
    t = parsed_args["t"]
    dt = parsed_args["dt"]
    plott = parsed_args["plot"]
    masses = parsed_args["masses"]
    names = parsed_args["names"]
    acc_func = parsed_args["acc_func"]
    int_func = parsed_args["int_func"]
    plot_dim = parsed_args["plot_dim"]
    return data, writefile, t, dt, plott, masses, names, acc_func, int_func, plot_dim
end

function dataSorter(data) #does black magic. Endrer fra defaultformatet i "npz", til det vi bruker i integrator
    items = Int64(length(data[1, 1, :]))
    dims = Int64(length(data[1, :, 1]))
    pos0 = zeros((items, dims))
    vel0 = zeros((items, dims))
    for i = 1:items
        pos0[i, :] = data[1, :, i]
        vel0[i, :] = data[2, :, i]
    end
    return items, vel0, pos0, dims
end
function forward_euler(vel0, pos0, t, dt, func, par, endvalue) 
"""
vel0 and pos0 should be (x, y, z) arrays with initial values for vel an pos
endvalues is a bool, return only resulting value after time t
otherwise return entire pos and vel array
"""
    
    #creating initial values
    len = Int64(round(t/dt))
    pl_len = length(vel0[:, 1])
    dims = length(vel0[1, :])
    vel = zeros((pl_len, dims, len))
    vel[:, :, 1] = vel0
    pos = zeros((pl_len, dims, len))
    pos[:, :, 1] = pos0
    ai = func(pos[:, :, 1], par)
    #integration loop
    for i = 2:len
        a = func(pos[:, :, i-1], par)
        vel[:, :, i] = vel[:, :, i-1] + a*dt #new vel
        pos[:, :, i] = pos[:, :, i-1] + dt*vel[:, :, i] #new pos
    end
    #return related stuff
    if endvalue
        return pos[:, :, end], vel[:, :, end]
    end
    return pos, vel
end

function velocity_verlet(vel0, pos0, t, dt, func, mas, endvalue)
"""
vel0 and pos0 should be [[planets],[x, y, z]] arrays with initial values for vel an pos
endvalues is a bool, return only resulting value after time t
otherwise return entire pos and vel array
"""

    #initial values
    len = Int64(round(t/dt))
    pl_len = length(vel0[:, 1])
    dims = length(vel0[1, :])
    vel = zeros((pl_len, dims, len))
    vel[:, :, 1] = vel0
    pos = zeros((pl_len, dims, len))
    pos[:, :, 1] = pos0
    ai = func(pos[:, :, 1], mas)
    #integration loop
    for i = 2:len
        aip = ai   #current acceleration
        pos[:, :, i] = pos[:, :, i-1] + dt*vel[:, :, i-1] + ((dt^2)/2)*aip  #new position
        ai = func(pos[:, :,  i], mas) #new acceleration based on new radius
        vel[:, :, i] = vel[:, :, i-1] + dt*(ai + aip)/2   #new velocity based on new acceleration
    end
    #return related stuff
    if endvalue
        return pos[:, :, end], vel[:, :, end]
    end
    return pos, vel
end

#6191060/160

function filewriter(dataArray, filename = "../prosjekt3/orbits.txt")
    #dataArray should be a 3d array, and is written to filename as a[:, 1]\n, a[:, 2]\n...
    f = open(filename,"w")
    for i = 1:Int64(length(dataArray)/3)
        npzwrite(f, dataArray[:, i])
    end
    close(f)
end


function integrate(vel0, pos0, t, dt, func, masses, int_func, endvalue = false)
    return int_func(vel0, pos0, t, dt, func, masses, endvalue)
end

function arrayify(thingy)
    thingy = split(thingy, r"'|[|]|,| ")
    filter!(e->e≠"",thingy)
    filter!(e->e≠"[",thingy)
    filter!(e->e≠"]",thingy)
    return thingy
end

function addMasses(namelist)
    #for some reason it was stupidly hard to read the masses from a file, si they are now hardcoded in
    maslib = Dict("SUN" => 1.989e30, "MERCURY" => 3.302e23, "VENUS" => 48.685e23, "EARTH" => 5.97219e24, "MARS" => 6.4171e23, "JUPITER" => 1.89813e27, "SATURN" => 5.6834e26, "URANUS" => 86.813e24, "NEPTUNE" => 102.413e24, "PLUTO" => 1.307e22)
    nem = arrayify(namelist)
    masses = zeros(Int64(length(nem)))
    for i=1:Int64(length(nem))
        masses[i] = maslib[nem[i]]
    end
    return masses/1.989e30
end

data, writefile, t, dt, plott, masses, names, acc_func, int_func, plot_dim= parse()
items, vel0, pos0, dims = dataSorter(data)
masses = addMasses(names)
acc_funcs = [aFunk, moreBodyFunc, moreBodyFunc_c]
acc_func = acc_funcs[acc_func]
int_funcs = [velocity_verlet, forward_euler]
int_func = int_funcs[int_func]
poss, vels = integrate(365.2242*vel0, pos0, t, dt, acc_func, masses, int_func)
plottify(poss, plot_dim, names)