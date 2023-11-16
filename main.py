import numpy as np
import matplotlib.pyplot as plt

# Define constants
dP_dx = -8
rho = 1
mu = 1
h = 1
points_y = 101
delta_t = 0.001
delta_y = 0.01
time = 0.025
time_points = int(time / delta_t)

# Create arrays
u = np.zeros((points_y, time_points))

# Define initial conditions
u[:,:] = 0

#Cranck-Nicolson method
for j in range(0, time_points-1):
    for n in range(0, points_y-1):
        u[n+1, j] = u[n,j] + delta_t * (-1/rho) * dP_dx + \
                    mu/2 * (((u[n+1,j+1] - 2*u[n+1,j] + u[n+1, j-1])/(delta_y**2)) + ((u[n,j+1] -2*u[n,j] + u[n,j-1])/(delta_y**2)))

plt.plot(u[:, 0], np.arange(1, points_y + 1))