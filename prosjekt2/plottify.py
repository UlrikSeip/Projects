import julia
import matplotlib.pytplot as plt
import scipy as sc

j = julia.Julia()
j.include(rotator.jl)

n = []
counter = [] 
a = sc.ones((i, i))
for i in range(100):
    newa, r, n_, counter_ = j.eval(rotate(a))
    n.append(n_)
    counter.append(counter_)

plt.plot(n, counter)
plt.show()
