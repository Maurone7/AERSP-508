import numpy as np
import matplotlib.pyplot as plt

# Constants
rho = 1.0
nu = 1.0
dPdx = -8.0
h = 1.0
Lx = 1.0
Ly = 1.0

# Discretization
Nx = 101  # Number of grid points in the x-direction
Ny = 101  # Number of grid points in the y-direction
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)
dt = 0.001  # Time step
T = 0.25  # Total simulation time

# Initialize velocity field
u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))

# Time-stepping loop
num_steps = int(T / dt)

for step in range(num_steps):
    # Crank-Nicolson scheme for u
    u_new = u.copy()
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            u_new[i, j] = u[i, j] + 0.5 * dt * (
                (nu / dx**2) * (u[i+1, j] - 2*u[i, j] + u[i-1, j] + u[i, j+1] - 2*u[i, j] + u[i, j-1])
                - (1 / rho) * (dPdx))
    u = u_new

# Plot the velocity profile
y = np.linspace(0, Ly, Ny)
plt.plot(u[int(Nx/2), :], y)
plt.xlabel('Velocity (u)')
plt.ylabel('y')
plt.title('Poiseuille Flow in the x-direction')
plt.grid()
plt.show()
