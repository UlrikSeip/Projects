import PyPlot
const plt = PyPlot

thingy = [1,2,3,4,6,7,6,4,3,6,3,7]
thingy2 = [1,2,3,3,4,6,7,7,8,4,3,7]

plt.plot3D(thingy, thingy2)
plt.show()