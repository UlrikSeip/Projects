function RungeKutta4(dt, f, y, t, var=[])
    """
    A function that implements one step of the Runge Kutta 4th order method.
    Takes the inputs for timestep dt, function f, current step y, current
    time t and a list var of other variables.
    Returns the next step yn.
    """
    dt2 = dt/2 #half the timestep
    k1 = f(t,y,var)
    k2 = f(t+dt2, y+k1*dt2, var)
    k3 = f(t+dt2, y+k2*dt2, var)
    k4 = f(t+dt, y+k3*dt, var)

    yn = y + (k1+k2+k3+k4)*dt/6
    return yn
end 
