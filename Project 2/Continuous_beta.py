import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import time

start_time = time.time()
fig, ax = plt.subplots()
falkner_skan_differential_equation = lambda eta, f, beta: [f[1], f[2], -f[0] * f[2] - beta * (1 - f[1] ** 2)]
secant = lambda x0, x1, f_1, f_0: x1 - (1 - f_1) * (x1 - x0) / ((1 - f_1) - (1 - f_0))
solve_falkner_skan_eqaution = lambda beta, guess: solve_ivp(falkner_skan_differential_equation, t_span=(0, 100), y0=[0, 0, guess], args=(beta,), t_eval=eta, method='BDF').y

min_beta, max_beta, N = -0.1988, 0.65, 40
beta_list, eta = list(np.linspace(min_beta, max_beta, N)), np.linspace(0, 100, int(100 / 0.001))
cmap = plt.get_cmap('nipy_spectral', N)
for beta in beta_list:
    beta = np.round(beta, 4)
    guess_1, guess_2 = 1, 0.1
    f_0 = solve_falkner_skan_eqaution(beta, guess_1)
    while abs(guess_2 - guess_1) > 1E-8:
        f_1 = solve_falkner_skan_eqaution(beta, guess_2)
        guess_1, guess_2, f_0 = guess_2, secant(guess_1, guess_2, f_1[1][-1], f_0[1][-1]), f_1
    print('Done', beta)
    f, f_prime, f_double_prime = f_1[0], f_1[1], f_1[2]
    plt.plot(f_prime[:5001], eta[:5001], label=r"$\beta={}$".format(beta), c=cmap(beta))

plt.xlabel(r"$\frac{u}{U_e}$"), plt.xlim(0, 1), plt.ylim(0, 5), plt.ylabel('$\eta$')
import matplotlib as mpl

norm = mpl.colors.Normalize(vmin=min_beta, vmax=max_beta)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, boundaries=np.arange(min_beta, max_beta, 0.05), ax=ax)
cbar.ax.set_title(r"$\beta$"), plt.show()

print("--- %s seconds --- Faster than MatLab" % (time.time() - start_time))
