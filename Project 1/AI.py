import numpy as np

# Parameters
rho = 1.0  # Fluid density
nu = 1.0  # Fluid viscosity
dPdx = -8.0  # Pressure gradient
h = 1.0  # Channel height

# Discretization parameters
Nx = 101  # Number of grid points in x direction
Ny = 101  # Number of grid points in y direction
dt = 0.001  # Time step

# Initialize velocity field
u = np.zeros((Ny, Nx))

# Boundary conditions
u[:, 0] = 0.0  # No-slip at bottom wall
u[:, Nx - 1] = 0.0  # No-slip at top wall
times_list = [0.025, 0.05, 0.1, 0.4]

# Crank-Nicolson scheme
for i in range(0, int(0.4 / dt) + 1):

  # Calculate intermediate velocities
  u_mid = u + dt / 2.0 * (-dPdx / rho + nu * np.gradient(u, axis=1))

  # Solve for new velocities
  A = np.eye(Ny) + dt * nu * (2 * np.diag(np.ones(Ny)) - np.diag(np.ones(Ny - 1), 1) - np.diag(np.ones(Ny - 1), -1))
  b = u_mid + dt * nu * np.gradient(u_mid, axis=1)
  u = np.linalg.solve(A, b)

# Plot velocity profile
import matplotlib.pyplot as plt

plt.plot(u[:, Ny // 2], np.arange(1, Nx + 1))
plt.xlabel('X')
plt.ylabel('Velocity')
plt.title('Velocity Profile of Unsteady Poiseuille Flow')
plt.show()
