using Random

function StI(N, a, S, I, dt)
    """
    Moves a number of people from S to I. This corresponds to a single time step.
    """
    counter = 0
    prob = a*S*I*dt/N
    rnum = rand()
    if rnum < prob
        counter += 1
    end
    return counter
end

function ItR(b, I, dt)
    """
    Moves a number of people from I to R. This corresponds to a single time step.
    """
    counter = 0
    prob = b*I*dt
    rnum = rand()
    if rnum < prob
        counter += 1
    end
    return counter
end

function RtS(c, R, dt)
    """
    Moves a number of people from R to S. This corresponds to a single time step.
    """
    counter = 0
    prob = c*R*dt
    rnum = rand()
    if rnum < prob
        counter += 1
    end
    return counter
end

function mcStep(a, b, c, S, I, R, N, dt)
    """
    Performs a single time step for all variables.
    """
    si = StI(N, a, S, I, dt)
    ir = ItR(b,I,dt)
    rs = RtS(c,R,dt)
    newS = S - si + rs
    newI = I - ir + si
    newR = R - rs + ir
    return newS, newI, newR
end