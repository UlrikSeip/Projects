function matrixmakinator(dim1, dim2, M, flip = false)
    #=
        Turns an array of arrays into a matrix.
    =#
    if flip
        newMatrix = zeros(dim2, dim1)
        for i = 1:dim1
            @inbounds for j = 1:dim2
                newMatrix[j, i] = M[i][j]
            end
        end
    else
        newMatrix = zeros(dim1, dim2)
        for i = 1:dim1
            @inbounds for j = 1:dim2
                newMatrix[i, j] = M[i][j]
            end
        end
    end
    return newMatrix
end

function statisticsinator(arrayarray)
    #=
        Does black magic, and finds average array along dim2 
        with the corresponding standard diviation at the final point
    =#
    dim1 = length(arrayarray)
    dim2 = length(arrayarray[1])
    M = matrixmakinator(dim1, dim2, arrayarray)
    Mavg = zeros(dim2)
    for i = 1:dim2
        Mavg[i] = sum(M[:, i])/dim1
    end
    stdDivM = stdm(M[:, end], Mavg[end])
    return Mavg, stdDivM
end