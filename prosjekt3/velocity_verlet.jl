include("integrator.jl")
include("afunks.jl")
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
        "acc_func"
            help = "a number representing which acceleration function you want to use; 1=aFunk, 2=acc_fs, 3=acc_nfs, 4=moreBodyFunc"
            required = false
            default = 1
            arg_type = Int64
    end
    return parse_args(args)
end

function Parse() #reads and returns initial info from file 
    parsed_args = parse_commandline()
    data = npzread(parsed_args["readfile"])
    writefile = parsed_args["writefile"]
    t = parsed_args["t"]
    dt = parsed_args["dt"]
    plott = parsed_args["plot"]
    masses = parsed_args["masses"]
    names = parsed_args["names"]
    acc_func = parsed_args["acc_func"]
    return data, writefile, t, dt, plott, masses, names, acc_func
end

function dataSorter(data) #does black magic. Endrer fra defaultformatet i "npz", til det vi bruker i integrator
    items = Int64(length(data[1, :, 1]))
    dims = Int64(length(data[1, 1, :]))
    pos0 = zeros((items, dims))
    vel0 = zeros((items, dims))
    print("items ")
    println(items)
    for i = 1:items
        pos0[i, :] = data[1, i, :]
        vel0[i, :] = data[2, i, :]
    end
    return items, vel0, pos0, dims
end

function plottify(items, names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptun", "Pluto"])
    counter = 1
    labels = split(names, r"'|[|]|,| ")
    filter!(e->e≠"",labels)
    filter!(e->e≠"[",labels)
    filter!(e->e≠"]",labels)
    println(labels)
    for i = 1:Int64(items/3)
        #plt.plot3D(poss[counter, :], poss[counter+1, :], poss[counter+2, :], label = labels[i])
        plt.plot(poss[counter, :], poss[counter+1, :], label = labels[i])
        counter += 3
    end
    plt.xlabel("pos x [AU]")
    plt.ylabel("pos y [AU]")
    #plt.zlabel("pos z [AU]")
    plt.axis("equal")
    plt.legend()
    plt.show()
end

acc_funcs = [aFunk, acc_fs, acc_nfs, moreBodyFunc]

data, writefile, t, dt, plott, masses, names, acc_func = Parse()    #creates variables for arguments
items, vel0, pos0, dims= dataSorter(data)   #sorts the data
mas = split(masses, r"'|\[|\]|,| ")
println(mas)
filter!(e->e≠"",mas)
filter!(e->e≠"[",mas)
filter!(e->e≠"]",mas)
println(mas)
#mas = [5.97219e24/2e30, 1.89813e27/2e30, 1]
poss, vels = velocity_verlet(365.2242*vel0, pos0, t, dt, acc_funcs[Int64(acc_func)], mas) #integrates
plottify(items, names)
#filewriter(poss, writefile) #writes to file