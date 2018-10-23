include("integrator.jl")
include("b_test.jl")
include("f_test.jl")
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
    return data, writefile, t, dt, plott, masses, names
end

function dataSorter(data) #does black magic. Endrer fra defaultformatet i "npz", til det vi bruker i integrator
    items = Int64(length(data[1, :, 1]))
    dims = Int64(length(data[1, 1, :]))
    pos0 = zeros((items, dims))
    vel0 = zeros((items, dims))
    for i = 1:items
        pos0[i, :] = data[1, i, :]
        vel0[i, :] = data[2, i, :]
    end
    return items, vel0, pos0, dims
end

function plottify(items)
    counter = 1
    labels = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptun", "Pluto"]
    for i = 1:Int64(items/3)
        """
        println(poss[counter, 1], " ", poss[counter+1, 1]," ", poss[counter+2, 1])
        println(vels[counter, 1], " ", vels[counter+1, 1]," ", vels[counter+2, 1])
        println(poss[counter, 2], " ", poss[counter+1, 2]," ", poss[counter+2, 2])
        println(vels[counter, 2], " ", vels[counter+1, 2]," ", vels[counter+2, 2])
        println()
        """
        plt.plot3D(poss[counter, :], poss[counter+1, :], poss[counter+2, :], label = labels[i])
        counter += 3
    end
    plt.xlabel("pos x [AU]")
    plt.ylabel("pos y [AU]")
    plt.zlabel("pos z [AU]")
    plt.axis("square")
    plt.legend()
    plt.show()
end

data, writefile, t, dt, plott, masses, names = parse()    #creates variables for arguments
items, vel0, pos0, dims= dataSorter(data)   #sorts the data
poss, vels = velocity_verlet(365.2242*vel0, pos0, t, dt, aFunk, []) #integrates
plottify(items)
#filewriter(poss, writefile) #writes to file