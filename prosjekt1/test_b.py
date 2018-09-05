import seaborn
from numpy import zeros, exp, linspace, ones
from matplotlib.pyplot import plot, show, title, xlabel, ylabel, legend, savefig

def rref(n):
    a = ones(n-1)
    a *= -1
    b = ones(n)
    b *= 2
    c = ones(n-1)
    c *= -1
    
    B = zeros(n)
    #v = zeros(n)
    b_ = zeros(n)
    k = zeros(n)
    B_ = zeros(n)
    x = linspace(0, 1, n)
    
    h = 1/(n+1)
    
    def f(x):
        return (100*exp(-10*x))
    
    u = 1 - (1 - exp(-10))*x - exp(-10*x)
    
    for i in range(n):
        B[i] = h**2*f(x[i])
    
        
    b_[0] = b[0]
    B_[0] = B[0]
    for i in range(1, n):
        temp = a[i-1]/b_[i-1]
        b_[i] = b[i] - temp*c[i-1]
        B_[i] = B[i] - temp*B_[i-1]
    
    k[n-1] = B_[n-1]
    i = n-2
    while i >= 0:
        k[i] = B_[i] - c[i]*k[i+1]/b_[i]
        i-=1
    return x, k, b_, u

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