using Random

function StI(N, a, S, I, dt)
    """
    Moves a number of people from S to I. This corresponds to a single time step.
    """
    prob = a*S*I*dt/N
    rnum = rand()
    if rnum < prob
        return 1
    end
    return 0
end

function ItR(b, I, dt)
    """
    Moves a number of people from I to R. This corresponds to a single time step.
    """
    prob = b*I*dt
    rnum = rand()
    if rnum < prob
        return 1
    end
    return 0
end

function RtS(c, R, dt)
    """
    Moves a number of people from R to S. This corresponds to a single time step.
    """
    prob = c*R*dt
    rnum = rand()
    if rnum < prob
        return 1
    end
    return 0
end

function nDeath(d, var, dt)
    prob = d*var*dt
    rnum = rand()
    if rnum < prob
        return 1
    end
    return 0
end

function iDeath(di, I, dt)
    prob = di*I*dt
    rnum = rand()
    if rnum < prob
        return 1
    end
    return 0
end

function vaccinated(f, S, dt)
    prob = f*S*dt
    rnum = rand()
    if rnum < prob
        return 1
    end
    return 0
end

function born(bi, S, dt)
    prob = bi*S*dt
    rnum = rand()
    if rnum < prob
        return 1
    end
    return 0
end

function mcStepB(a, b, c, S, I, R, N, dt)
    """
    Performs a single time step for all variables, with the simple version of the model.
    """
    si = StI(N, a, S, I, dt)
    ir = ItR(b,I,dt)
    rs = RtS(c,R,dt)
    newS = S - si + rs
    newI = I - ir + si
    newR = R - rs + ir
    return newS, newI, newR
end

function mcStep(A, omega, t, a0, b, c, d, di, bi, f, S, I, R, N, dt)
    a = A*cos(omega*t) + a0
    """
    Performs a single time step for all variables.
    """
    si = StI(N, a, S, I, dt)
    ir = ItR(b,I,dt)
    rs = RtS(c,R,dt)

    dS = nDeath(d, S, dt)
    dI = nDeath(d, I, dt)
    dR = nDeath(d, R, dt)

    diI = iDeath(di, I, dt)
    vac = vaccinated(f, S, dt)
    bS = born(bi, S, dt)

    newS = S - si + rs - dS - vac + bS
    newI = I - ir + si - dI - diI
    newR = R - rs + ir - dR + vac
    return newS, newI, newR
end
