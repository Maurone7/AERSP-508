import numpy as np

B_n = lambda n: 8 * ((np.pi * n)*np.sin(np.pi*n) + 2 * np.cos(np.pi * n) - 2) / (np.pi ** 3 * n ** 3)
n_plot = []
B_n_plot = []

for n in range(1, 500):
    if n % 2 != 0:
        n_plot.append(n)
        B_n_plot.append(B_n(n))
        if abs(B_n(n))< 1e-8:
            n_max = n
            print(n_max)
            break
