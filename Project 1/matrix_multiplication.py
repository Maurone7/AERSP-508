import numpy as np

matrix_1 = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
matrix_2 = np.array([[1], [2], [3]])

print(np.shape(matrix_2))

print(matrix_1 * matrix_2)
print(np.linalg.solve(matrix_1, matrix_2))