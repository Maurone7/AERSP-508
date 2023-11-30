import numpy as np
import matplotlib.pyplot as plt
# Define falker-skan equation
def falker_skan_equation(beta, f, f_prime, f_double_prime, f_triple_prime):
    return f_triple_prime + f * f_double_prime + beta * (1 - f_prime**2)

# Values of eta
eta = np.linspace(0, 100, int(100/0.001))

# make a list of m values from 0 to 1 as a list
m_list = np.linspace(0, 1, 11)
beta = lambda m: 2*m/(1+m)
beta_list = [beta(m) for m in m_list]

def rk4(fn, y0, t, dt):
    """
    RK4 integrator

    Parameters:
    - fn: the function to be integrated
    - y0: the initial value of the state vector
    - t: the initial time
    - dt: the time step

    Returns:
    - y1: the state vector at the next time step
    """

    k1 = fn(t, y0)
    k2 = fn(t + 0.5 * dt, y0 + 0.5 * dt * k1)
    k3 = fn(t + 0.5 * dt, y0 + 0.5 * dt * k2)
    k4 = fn(t + dt, y0 + dt * k3)

    y1 = y0 + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return y1


