function dS(t,S,var)
    """
    A function that finds the change in susceptible individuals in the SIRS
    model.
    """
    c,R,a,I,N = var
    s = c*R - a*S*I/N
    return s
end

function dI(t,I,var)
    """
    A function that finds the change in infected individuals in the SIRS
    model.
    """
    a,S,N,b = var
    i = - b*I + a*S*I/N
    return i
end

function dR(t,R,var)
    """
    A function that finds the change in recovered individuals in the SIRS
    model.
    """
    c,I,b = var
    r = b*I - c*R
    return r
end

function dSvd(t,S,var)
    """
    A function that finds the change in susceptible individuals in the SIRS
    model, when acounting for vital dynamics.
    """
    c,R,a,I,N,d,e = var
    s = c*R - a*S*I/N - d*S + e*N
    return s
end

function dIvd(t,I,var)
    """
    A function that finds the change in infected individuals in the SIRS
    model, when acounting for vital dynamics.
    """
    a,S,N,b,d,di = var
    i =  a*S*I/N - b*I - d*I - di*I
    return i
end

function dRvd(t,R,var)
    """
    A function that finds the change in recovered individuals in the SIRS
    model, when acounting for vital dynamics.
    """
    c,I,b,d = var
    r = b*I - c*R - d*R
    return r
end
