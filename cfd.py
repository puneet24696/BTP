import numpy as np
import matplotlib.pyplot as plt

rho = 14.33 * 9.81
T = 10**6
f = 30

min_x = 0
max_x = 200
min_t = 0
max_t = 5
N = 200
dx = (max_x - min_x) / (N - 1)
CFL = 0.5
dt = CFL * rho * dx * dx / T

x, t = np.mgrid[min_x:max_x:dx, min_t:max_t:dt]
u = x.copy()
print u, x

u[:,:] = 0
u[0,:] = 0
u[:,0] = 0
u[:,1] = 0
u[-1,:] = 0

for j in range(2, len(u[0, :])):
    for i in range(1, len(u[:,0]) - 1):
        u[i,j] = (T * dt * dt / rho * dx * dx) * (u[i+1, j-2] - 2 * u[i, j-2] + u[i-1, j-2]) + (f * dt * dt / rho) +  2 * u[i, j-1] - u[i, j-2]

print u

for i in range(0,len(u[0,:]),1000):
    plt.plot(x[:,0], u[:,i])
    plt.ylim([0, 1])
    plt.title("Vibrating String at time "+str(i*5/70000.)+"s")
    plt.xlabel("length of sting(m)")
    plt.ylabel("Deflection(m)")
    plt.savefig("fig:" + str(i))
    plt.waitforbuttonpress(0.00001)
    plt.cla()

