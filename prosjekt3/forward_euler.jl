include("integrator.jl")
include("b_test.jl")
import PyPlot
const plt = PyPlot
using ArgParse
using DelimitedFiles
using NPZ
using PyPlot


"""
Simply run the file with the required command line arguments to integrate the supplied values.
The values should be stored in a file with the format:

####
9.413801075750535E-01 3.379019986046322E-01 -9.334104672733438E-05      #initial positions of object 1 x, y z
-5.994522787486753E-03 1.617377250092178E-02 -1.732657683299539E-07     #initial velocity of object 1 x, y z
1.375357774690336E+00 -1.627517936385902E-01 -3.738675132962930E-02     #initial positions of object 2 x, y z
2.243217631343320E-03 1.508628660091320E-02 2.610262676274213E-04       #initial velocity of object 2 x, y z
####

I am still not sure how one should go about supplying a custom force function
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
    return data, writefile, t, dt, plott
end

function dataSorter(data) #does black magic. Endrer fra default formatet til "npz", til det vi bruker i integrator
    items = Int64(length(data[1, :, 1]))
    dims = Int64(length(data[1, 1, :]))
    pos0 = zeros((items, dims))
    vel0 = zeros((items, dims))
    for i = 1:items
        for j = 1:dims
            pos0[i, j] = data[1, i, j]
            vel0[i, j] = data[2, i, j]
        end
    end
    return items, vel0, pos0
end

function plottify(items)
    counter = 1
    for i = 1:items
        println(poss[counter], poss[counter+1], poss[counter+2])
        plt.plot3D(poss[counter], poss[counter+1], poss[counter+2])
        counter += 3
    end
    plt.show()
end

data, writefile, t, dt, plott= parse()    #creates variables for arguments
items, vel0, pos0 = dataSorter(data)   #sorts the data
poss, vels = forward_euler(365.2422*vel0, pos0, t, dt, aFunk, []) #integrates
println(size(poss))
if plott
    plottify(items)
end
filewriter(poss, writefile) #writes to file