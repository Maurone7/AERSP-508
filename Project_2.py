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

for betas in beta_list:
    