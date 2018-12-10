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

function dSda(t,S,var)
    """
    A function that finds the change in susceptible individuals in the SIRS
    model when a can vary.
    """
    c,R,A,a0,om,I,N = var
    a = da(A,a0,om,t)
    s = c*R - a*S*I/N
    return s
end

function dIda(t,I,var)
    """
    A function that finds the change in infected individuals in the SIRS
    model when a can vary.
    """
    A,a0,om,S,N,b = var
    a = da(A,a0,om,t)
    i = - b*I + a*S*I/N
    return i
end

function da(A,a0,om,t)
    """
    A function that modeles a as a sin-function. 
    """
    #we assume we start with average transmission rate
    return A*sin.(t*om) + a0
    #return A*cos.(t*om) + a0
end

function dSf(t,S,var)
    """
    A function that finds the change in susceptible individuals in the SIRS
    model, when accounting for vaccination, represented by f.
    """
    c,R,a,I,N,f_funk,fvar = var
    f = f_funk(t,fvar)
    s = c*R - a*S*I/N - f*S
    return s
end

function dRf(t,R,var)
    """
    A function that finds the change in recovered individuals in the SIRS
    model, when accounting for vaccination, represented by f.
    """
    c,I,b,S,f_funk,fvar = var
    f = f_funk(t,fvar)
    r = b*I - c*R + f*S
    return r
end

function df(t,fvar)
    """
    A function that modeles f as a linear function.
    """
    f0,Df = fvar
    return f0 + t*Df
end
