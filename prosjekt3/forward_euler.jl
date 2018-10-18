include("integrator.jl")
include("b_test.jl")
using ArgParse
using DelimitedFiles

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
    end
    return parse_args(args)
end

function parse() #reads and returns initial info from file 
    parsed_args = parse_commandline()
    data = readdlm(parsed_args["readfile"])
    writefile = parsed_args["writefile"]
    t = parsed_args["t"]
    dt = parsed_args["dt"]
    return data, writefile, t, dt
end

function dataSorter(data) #does black magic. Endrer fra default formatet til "readdlm", til det vi bruker i integrator
    items = Int64(length(data[:, 1])/2)
    varsPerItem = Int64(length(data[1, :]))
    pos0 = zeros((items, varsPerItem))
    vel0 = zeros((items, varsPerItem))
    for i = 1:items
        pos0[i, :] = data[i, :]
        vel0[i, :] = data[i*2, :]
    end
    return vel0, pos0
end

data, writefile, t, dt = parse()    #creates variables for arguments
vel0, pos0 = dataSorter(data)   #sorts the data
poss, vels = forward_euler(365.2422*vel0, pos0, t, dt, aFunk) #integrates
filewriter(poss, writefile) #writes to file