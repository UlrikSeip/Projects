import seaborn
from numpy import zeros, exp, linspace, ones
from matplotlib.pyplot import plot, show, title, xlabel, ylabel, legend, savefig

def plot_rref(n):
    x, k, b_, u = rref(n) 
    
    plot(x, k)
    leg.append("n: " + str(n))
    return x, u

leg = []
plot_rref(10)
plot_rref(100)
x, u = plot_rref(1000)
plot_rref(10000)


plot(x, u)
xlabel("x"); ylabel("f(x)")
leg.append("exact")
legend(leg)
title("Plot of the aprox. and exact solutions")
savefig("test_b.pdf")
show()

def plot_rref(n):
    x, k, b_, u = rref(n) 
    
    plot(x, k)
    leg.append("n: " + str(n))
    return x, u

leg = []
plot_rref(10)
plot_rref(100)
x, u = plot_rref(1000)

plot(x, u)
xlabel("x"); ylabel("f(x)")
leg.append("exact")
legend(leg)
title("Plot of the aprox. and exact solutions")
savefig("test_b.pdf")
show()