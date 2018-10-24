using LinearAlgebra

function aFunk(pos, mas)
    
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

function twoBodyFunc(pos1, pos2, mas2)
    #G = 6.67408e-11
    G = -4*(pi^2)
    #auToM = 149.60*10^9
    items = Int64(length(pos1[:, 1, 1]))
    dims = Int64(length(pos1[1, :, 1]))
    r = norm(pos2-pos1)
    th = (pos2-pos1)/r
    a = G*mas2*th/((r)^2)
    return a
end

function moreBodyFunc(pos, mas) #positions and masses in arrays
    items = Int64(length(pos[:, 1, 1]))
    dims = Int64(length(pos[1, :, 1]))
    a = zeros((items, dims))
    for i = 1:items                     #loops over and finds a on i from j
        planpos = zeros((1, dims, 1))
        planpos[1, :] = pos[i, :]
        aa = zeros(dims)
        aa += aFunk(planpos, mas[i])[1, :] #a from sun
        for j = 1:items
            if i != j
                aa += twoBodyFunc(pos[i, :], pos[j, :], mas[i]) #a from planet j
            end
        end
        a[i, :] = aa
    end
    return a
end


function moreBodyFunc_c(pos, mass)
end