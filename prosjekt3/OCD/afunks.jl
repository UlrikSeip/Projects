using LinearAlgebra

function aFunk(pos, mass)
    
    """
    Basic acceleration function with all the bodies orbiting the sun, and only the sun
    Thus could probalby be made prettier, but array operations are tricky in julia (at least i find them to be..)
    Returns:
        [float array[items, dims]] -- [acceleration in the dimensions of the given positions]
    """

    sunMass = 1
    #G = 6.67408e-11
    G = -4*(pi^2)
    #auToM = 149.60*10^9
    items = Int64(length(pos[:, 1, 1]))
    dims = Int64(length(pos[1, :, 1]))
    r = zeros(items)
    for i = 1:items
        r[i] = norm(pos[i, :])
    end
    th = zeros(items, dims)
    for i = 1:items
        th[i, :] = pos[i, :]/r[i]
    end
    a = zeros(items, dims)
    for i = 1:items
        a[i, :] = G*sunMass*th[i, :]/((r[i])^2)
    end
    return a
end

function twoBodyFunc(pos1, pos2, mass2)
    sunMass = 1
    #G = 6.67408e-11
    G = -4*(pi^2)
    #auToM = 149.60*10^9
    items = Int64(length(pos1[:, 1, 1]))
    dims = Int64(length(pos1[1, :, 1]))
    r = zeros(items)
    for i = 1:items
        r[i] = norm(pos[i, :])
    end
    th = zeros(items, dims)
    for i = 1:items
        th[i, :] = pos[i, :]/r[i]
    end
    a = zeros(items, dims)
    for i = 1:items
        a[i] = G*sunMass*th[i]/((r[i])^2)
    end
    return a
end

function moreBodyFunc(pos, mass) #positions and masses in arrays
    items = Int64(length(pos)/3)
    poss = zeros(3, items)
    for i=1:Int64(items)
        poss[1, i] = pos[i]
        poss[2, i] = pos[i+1]
        poss[3, i] = pos[i+2]
    end
    mass = [2.986095e-06, 9.490650e-04]
    append!(pos, [0,0,0])
    append!(mass, 1)
    a = zeros((3, items-1))
    for j = 1:items-1
        aa = 0
        for i = 1:items-1
            if i != j
                aa += twoBodyFunc(poss[:, j], poss[:, i], mass[i])
            end
        end
        a[j] = aa
    end
    return a
end


function moreBodyFunc_c(pos, mass)
end