using Printf
using Statistics
const stat = Statistics

const hbar = 4.135667662e-15                      #eVs
const omega = 1                                   #rad/s
const m = 0.5110e-6                               #eV
const k = m*omega^2
const alpha = (hbar^2/(m*k))^0.25
const l_factor = (2*m*alpha^2)/(hbar^2)

const rho_min = 1e-6
rho_max = 10

function E(n)
    #=
    Function to analytically calculate the
    energy levels of the Quantum System
    Takes arguments:
    n - energy state number (Int64)
    =#
    return hbar*omega*(2*n + 1.5)
end

function SQM_Tridiag(N)
    #=
    Function to create a Tridiagonal matrix with
    increasing values representing potential energy for a single electron
    along the diagonal.
    Takes the argument:
    N - Dimension of the NxN matrix (Int64)
    =#

    h = (rho_max - rho_min)/N
    M = zeros(Float64, (N,N))
    e = -1/(h^2)

    #Setting boundary cases
    M[1,1] = 2/(h^2) + rho_min^2
    M[1,2] = e
    M[N,N] = 2/(h^2) + rho_max^2
    M[N,N-1] = e

    for i=2:N-1
        rho = rho_min + i*h
        M[i,i] = 2/(h^2) + rho^2
        M[i,i-1] = e
        M[i,i+1] = e
    end
    return M
end

function TQM_Tridiag(N, wr)
    #=
    Function to create a Tridiagonal matrix with
    increasing values representing potential energy for two electrons
    along the diagonal.
    Takes the arguments:
    N - Dimension of the NxN matrix (Int64)
    wr - Is the scaled 'frequency' of the system (Float64)
    =#

    h = (rho_max - rho_min)/N
    M = SQM_Tridiag(N)

    #Setting boundary cases
    M[1,1] = 2/(h^2) + (wr*rho_min)^2 + 1/rho_min
    M[N,N] = 2/(h^2) + (wr*rho_max)^2 + 1/rho_max

    for i=2:N-1
        rho = rho_min + i*h
        M[i,i] = 2/(h^2) + (wr*rho)^2 + 1/rho

    end
    return M
end

function jacobi_method(A, eps=1e-5)
    #=
    Function for finding the eigenvalues of the symmetrical matrix A
    using Jacobi's rotation algorithm (rotation is implemented in the function rotate).
    Takes the arguments
    A - The matrix on the form [(NxN)] with elements with dtype Float64
    eps - Epsilon, the limit for when the program will consider the matrix element to be zero (Float64)
    =#
    N = size(A)[1]
    R = zeros((N,N))

    for i=1:N
        for j=1:N
            i == j ? R[i,j] = 1.0 : R[i,j] = 0.0
        end
    end

    max_number_iterations = N^3
    iterations = 0
    max_offdiag, k, l = maxoffdiag(A, 1,  1, N)

    while (abs(max_offdiag) > eps) && (iterations < max_number_iterations)
        max_offdiag, k, l = maxoffdiag(A, k, l, N)
        rotate(A, R, k, l, N)
        iterations += 1
    end
    return A
end

function maxoffdiag(A, k, l, N)
    #=
    A simple function for finding the maximum value of the matrix A
    which isn't placed on the diagonal:
    A - a symmetrical NxN matrix
    N - dimension of A
    k,l - previous maximum indices, in case none are found
    =#
    max = 0.0
    for i=1:N
        for j=i+1:N
            m = abs(A[i,j])
            if (m > max)
                max = m
                k = j
                l = i
            end
        end
    end
    return max, k, l
end

function rotate(A, R, k, l, N)
    #=
    This is where the magic happens
    Implementation of the rotation in Jacobi's rotation algorithm
    takes the arguments
    A - Symmetrical matrix with dimensions NxN and dtype Float64
    R - An NxN identity matrix
    k,l - The x- and y-indexes of the largest non-diagonal element
    N - Dimension of A
    =#
    c = 1.0
    s = 0.0

    if (A[k,l] != 0.0)
        tau = (A[l,l] - A[k,k])/(2*A[k,l])
        if (tau > 0)
            t = -tau + sqrt(1.0 + tau^2)
        else
            t = -tau - sqrt(1.0 + tau^2)
        end
        c = 1/sqrt(1.0 + t^2)
        s = c*t
    end
    a_kk = A[k,k]
    a_ll = A[l,l]
    s2 = s^2
    c2 = c^2
    A[k,k] = a_kk*c2 - 2.0*c*s*A[k,l] + a_ll*s2
    A[l,l] = a_kk*s2 + 2.0*c*s*A[k,l] + a_ll*c2
    A[k,l] = 0.0
    A[l,k] = 0.0
    for i=1:N
        if (i != k) && (i != l)
            a_ik = A[i,k]
            a_il = A[i,l]
            A[i,k] = c*a_ik - s*a_il
            A[k,i] = A[i,k]
            A[i,l] = c*a_il + s*a_ik
            A[l,i] = A[i,l]
        end
        r_ik = R[i,k]
        r_il = R[i,l]
        R[i,k] = c*r_ik - s*r_il
        R[i,l] = c*r_il + s*r_ik
    end
end

function print_matrix(M)
    N = size(M)[1]
    for i=1:N
        for j=1:N
            print("$(@sprintf("%7.2f ", M[i,j]))")
        end
        println()
    end
end

function get_eigenvalues(A)
    #=
    A wrapper that applies the Jacobi method to an arbitrary symmetrical
    matrix A, and then extracts the diagonal and sorts it from low to high.
    =#

    #Getting the matrix dimension
    N = size(A)[1]

    #Applying the Jacobi algorithm
    D = jacobi_method(A)

    #Initializing the eigenvalue array
    eigenvalues = zeros(N)
    #Inserts the diagonals into the blank eigenvalue array
    for i=1:N
        eigenvalues[i] = D[i,i]
    end
    #Sorts the eigenvalues from low to high
    return sort(eigenvalues)
end

function get_analytical_eigenvalues(N)
    #=
    Calculates the analytical eigenvalues for N electron energy levels in
    atomic Hydrogen
    =#
    #Initializing the eigenvalue array
    eigenvalues = zeros(N)
    #Calculates each eigenvalue using the product of E_n and lambda
    for n=1:N
        eigenvalues[n] = l_factor*E(n-1)
    end
    #Sorts the eigenvalues from low to high
    return sort(eigenvalues)
end

#TEST FUNCTIONS

function test_TQM_Tridiag(tol = 1e-4)
    A = [10, 6, 5, 2]
    B = [1, 0.1, 0.5, 0.01]
    C = [
    [1.0e6 -1 0 0 0 0 0 0 0 0; -1 6.5 -1 0 0 0 0 0 0 0; 0 -1 11.3333 -1 0 0 0 0 0 0;
    0 0 -1 18.25 -1 0 0 0 0 0; 0 0 0 -1 27.2 -1 0 0 0 0; 0 0 0 0 -1 38.1667 -1 0 0 0;
    0 0 0 0 0 -1 51.1429 -1 0 0; 0 0 0 0 0 0 -1 66.125 -1 0; 0 0 0 0 0 0 0 -1 83.1111 -1;
    0 0 0 0 0 0 0 0 -1 102.1],
    [1.0e6 -0.36 0 0 0 0; -0.36 1.13111 -0.36 0 0 0; 0 -0.36 1.17 -0.36 0 0;
    0 0 -0.36 1.31444 -0.36 0; 0 0 0 -0.36 1.53444 -0.36; 0 0 0 0 -0.36 1.82],
    [1.0e6 -0.25 0 0 0; -0.25 4.75 -0.25 0 0; 0 -0.25 9.66667 -0.25 0;
    0 0 -0.25 16.625 -0.25; 0 0 0 -0.25 25.6],
    [1.0e6 -0.04; -0.04 0.19]
    ]
    for i=1:length(A)
        D = TQM_Tridiag(A[i], B[i])
        for j=1:A[i]
            for k=1:A[i]
                if abs((C[i][j,k] - D[j,k])/C[i][j,k]) > tol
                    error("Function \"SQM_Tridiag\" returns incorrect value for N = $(A[i]) and wr = $(B[i]) at:
                    \nIndex: ($j,$k)\nExpected: $(C[i][j,k]), got: $(D[j,k])\n")
                end
            end
        end
    end
end

function test_SQM_Tridiag(tol = 1e-4)
    A = [10, 6, 5, 2]
    B = [
    [2 -1 0 0 0 0 0 0 0 0; -1 6 -1 0 0 0 0 0 0 0; 0 -1 11 -1 0 0 0 0 0 0;
    0 0 -1 18 -1 0 0 0 0 0; 0 0 0 -1 27 -1 0 0 0 0; 0 0 0 0 -1 38 -1 0 0 0;
    0 0 0 0 0 -1 51 -1 0 0; 0 0 0 0 0 0 -1 66 -1 0; 0 0 0 0 0 0 0 -1 83 -1;
    0 0 0 0 0 0 0 0 -1 102],
    [0.72 -0.36 0 0 0 0; -0.36 11.8311 -0.36 0 0 0; 0 -0.36 25.72 -0.36 0 0;
    0 0 -0.36 45.1644 -0.36 0; 0 0 0 -0.36 70.1644 -0.36; 0 0 0 0 -0.36 100.72],
    [0.5 -0.25 0 0 0; -0.25 16.5 -0.25 0 0; 0 -0.25 36.5 -0.25 0;
    0 0 -0.25 64.5 -0.25; 0 0 0 -0.25 100.5],
    [0.08 -0.04; -0.04 100.08]
    ]
    for i=1:length(A)
        C = SQM_Tridiag(A[i])
        for j=1:A[i]
            for k=1:A[i]
                if abs(C[j,k] - B[i][j,k]) > tol
                    error("Function \"SQM_Tridiag\" returns incorrect value for N = $(A[i]) at:
                    \nIndex: ($j,$k)\nExpected: $(B[i][j,k]), got: $(C[j,k])\n")
                end
            end
        end
    end
end

function test_maxoffdiag(tol = 1e-4)
    A = [
    [2 1 0 0 0; 1 2 1 0 0; 0 1 2 1 0; 0 0 1 2 1; 0 0 0 1 2],
    [2 1; 1 2],
    [3 -6 0; -6 0 6; 0 6 -3],
    [9 2 5 2; 2 2 6 8; 5 6 -10 -1; 2 8 -1 6],
    [9 -9 -5 4 -7 -5; -9 2 -8 2 7 7; -5 -8 0 6 -2 4; 4 2 6 -2 4 9; -7 7 -2 4 4 0; -5 7 4 9 0 7],
    [-7 3 -2 7; 3 6 -3 7; -2 -3 -1 3; 7 7 3 7]
    ]

    test_max = [1, 1, 6, 8, 9, 7]
    test_l = [1, 1, 1, 2, 1, 1]
    test_k = [2, 2, 2, 4, 2, 4]

    for i=1:length(test_max)
        max, k, l = maxoffdiag(A[i], 0, 0, size(A[i])[1])
        if abs(max - test_max[i]) > tol
            error("\nfunction \"maxoffdiag\" does not return the correct maximum value
            \nExpected value: $(test_max[i]), got: $max.")
        elseif (k != test_k[i]) || (l != test_l[i])
            error("\nfunction \"maxoffdiag\" does not return the maximum indices
            \nExpected indices ($(test_k[i]), $(test_l[i])), got ($k, $l)\n")
        end
    end
end

function test_jacobi_method(tol = 1e-4)
    #=
    Checks that the difference between the true eigenvalues and the
    eigenvalues calculated in this program has a standard deviation within
    the given tolerance, for the given matrices A.  B are the true
    eigenvalues for the corresponding A matrices.
    =#

    #Test Matrices
    A = [
    [2 1 0 0 0; 1 2 1 0 0; 0 1 2 1 0; 0 0 1 2 1; 0 0 0 1 2],
    [1 2; 2 1],
    [3 -6 0; -6 0 6; 0 6 -3],
    [9 2 5 2; 2 2 6 8; 5 6 -10 -1; 2 8 -1 6],
    [9 -9 -5 4 -7 -5; -9 2 -8 2 7 7; -5 -8 0 6 -2 4; 4 2 6 -2 4 9; -7 7 -2 4 4 0; -5 7 4 9 0 7],
    [-7 3 -2 7; 3 6 -3 7; -2 -3 -1 3; 7 7 3 7]
    ]

    #True Eigenvalues
    B = [
    [0.2679, 1.0000, 2.0000, 3.0000, 3.7321],
    [-1, 3],
    [-9, 0, 9],
    [-14.1067, -1.7515,  7.8452, 15.0129],
    [-15.6343, -10.5784, 1.6880, 6.7285, 13.1854, 24.6108],
    [-10.9113, -3.4031, 3.4993, 15.8151]
    ]

    #Loops through each A[n]/B[n] and runs A[n] through our eigenvalue solver
    for n=1:size(A)[1]
        #Calculate the eigenvalues for the current matrix A[n]
        C = get_eigenvalues(convert(Array{Float64}, A[n]))

        #Standard deviation of the difference between C and B[n]
        std = stat.std(abs.(sort(B[n])-C))

        #Raises an error if the standard deviation is too high
        if std > tol
            error("\nFunction \"jacobi_method\" returns incorrect (within a tolerance of $(tol)) eigenvalues for matrix:
            \n$(A[n])\n\nExpected eigenvalues: $(B[n]), calculated eigenvalues: $C\n")
        end
    end
end

function test_get_analytical_eigenvalues(tol = 1e-4)
    A = [3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71,
    75, 79, 83, 87, 91, 95, 99, 103, 107, 111, 115, 119, 123, 127, 131, 135,
    139, 143, 147, 151, 155, 159, 163, 167, 171, 175, 179, 183, 187, 191, 195,
    199, 203, 207, 211, 215, 219, 223, 227, 231, 235, 239, 243, 247, 251, 255,
    259, 263, 267, 271, 275, 279, 283, 287, 291, 295, 299, 303, 307, 311, 315,
    319, 323, 327, 331, 335, 339, 343, 347, 351, 355, 359, 363, 367, 371, 375,
    379, 383, 387, 391, 395, 399.0]
    B = get_analytical_eigenvalues(100)
    for i=1:length(A)
        if abs(A[i] - B[i]) > tol
            error("Function \"get_analytical_eigenvalues\" failed at energy level n = $i
            \nExpected eigenvalue: $(A[i]), got: $(B[i])\n")
        end
    end
end

function run_tests(tol = 1e-4)
    test_SQM_Tridiag(tol)
    test_TQM_Tridiag(tol)
    test_maxoffdiag(tol)
    test_get_analytical_eigenvalues(tol)
    test_jacobi_method(tol)
end

#MAIN PROGRAM

function main(experiment)
    #=
    A function to run the different numerical experiments
    the variable experiment can be set to the following (string values):
    single - Finds the energy eigen values of a single electron in a harmonic oscillator potential
    pair - same as single, but with two electrons

    =#
    if experiment=="single"
        N = 400
        M = SQM_Tridiag(convert(Int64, N))
        @time eigenvalues = get_eigenvalues(M)
        analytical_eigenvalues = get_analytical_eigenvalues(N)
        diff = eigenvalues-analytical_eigenvalues
        l_mean = stat.mean(diff)
        l_std = stat.std(diff)
        for i=1:N
            println("$(@sprintf("%5.g", analytical_eigenvalues[i]))  $(@sprintf("%5.g", eigenvalues[i]))")
        end
    elseif experiment=="pair"
        omega = [0.01, 0.5, 1, 5]
        for wr in omega
            N = 400
            M = TQM_Tridiag(convert(Int64, N), wr)
            @time eigenvalues = get_eigenvalues(M)
            println(eigenvalues[1:7])
        end
    end
end

"""
experiment = "single"
run_tests()
main(experiment)
"""