function aFunk(pos, par)
    mass = sunM
    #G = 6.67408e-11
    G = -4*(pi^2)
    #auToM = 149.60*10^9
    r = norm(pos)
    th = pos/r
    a = G*mass*th/((r)^2)
    return a
end

function twoBodyFunc(pos1, pos2, mass2)
    mass = mass2
    #G = 6.67408e-11
    G = -4*(pi^2)
    #auToM = 149.60*10^9
    r = norm(pos2-pos1)
    th = (pos2-pos1)/r
    a = G*mass*th/((r)^2)
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
