include("rotator.jl")

function test_maxKnotL()
    test_matrix = ones(Float64, 10, 10)
    test_matrix[1, 4] = 3
    test_matrix[4, 1] = 3
    k, l = maxKnotL(test_matrix)
    if test_matrix[k, l] < 3
        error("maxKnotL returns wrong value")
    end
end

function test_eigenvalues()
    #creating testable values
    n = Int64(5)
    a = zeros(Float64, n, n)
    for i in range(1, step = 1, length = n)
        for j in range(1, step = 1, length = n)
            if ((j == i+1) || (j == i-1))
                a[i, j] = 2                     # b = 2
            end
            if (i == j)
                a[i, j] = 5                     # a = 5
            end
        end
    end
    test_values = [5-2*sqrt(3), 3, 5, 7, 5+2*sqrt(3)]
    #calculating values
    a_, thing1, thing2, thing3, thing4 = rotate(a, 1e-5)
    tol = 1e-3
    #sorting values into usable format
    values = zeros(5)
    for i = 1:5
        values[i] = a_[i,i]
    end
    values = sort(values)
    #tests each element
    for i = n:5
        if abs(values[i]-test_values[i]) > tol
            error("Eigentest failed: expected: $test_values, got $values")
        end
    end
end

test_eigenvalues()