import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

beta_list = [-0.1988, -0.1, -0.05, 0, 0.2, 0.5]
eta = np.linspace(0, 5, 5001)

def falkner_skan_differential_equation(eta, f, beta):
    return [f[1], f[2], -f[0] * f[2] - beta * (1 - f[1] ** 2)]

for beta in beta_list:
    f_init = [0, 0, 1]
    f, f_prime, f_double_prime = solve_ivp(falkner_skan_differential_equation, t_span=(0, 5), y0=f_init, args=(beta,), t_eval=eta, method='RK45').y
    while f_prime[-1] > 1:
        f_init[2] -= 0.001
        f, f_prime, f_double_prime = solve_ivp(falkner_skan_differential_equation, t_span=(0, 5), y0=f_init, args=(beta,), t_eval=eta, method='RK45').y
    plt.plot(f_prime[:len(eta)], eta, label=r"$\beta={}$".format(beta))
plt.xlabel('f'), plt.xlim(0,1), plt.ylim(0,5), plt.ylabel('$\eta$'), plt.legend(), plt.show()
