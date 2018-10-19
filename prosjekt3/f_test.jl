#include("integrator.jl")
#using PyCall
#using PyPlot
using LinearAlgebra
#using Math
#@pyimport matplotlib.pyplot as plotter
#@pyimport numpy as np
import LinearAlgebra: norm


function acc_fs(pos, s_mass)
    """
    Finds the acceleration for mulitple planets when the sun is fixed in the centre
    """
    G = -4*(pi^2)
    nr_pl = length(pos[:,:])/3
    
    pl = []
    i = 1
    while i < length(pos[:,:])
        append!(pl, pos[i:i+2, :])
        i+=3
    end
   
    println(pl)
    println(length(pl[1,:]))
    println(nr_pl)

    r = norm(pos)
    th = pos/r
    a = G*s_mass*th/((r)^2)
    return a
end

test = [[1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1]]
acc_fs(test, 1)

