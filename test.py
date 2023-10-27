import unittest
from functools import partial

import numpy as np
from typing import List


def lu_decomposition_determinant(matrix):
    def _lu_decomposition(matrix):
        n = len(matrix)
        lower = [[0.0] * n for _ in range(n)]
        upper = [[0.0] * n for _ in range(n)]

        # Decomposition
        for i in range(n):
            # Upper Triangular
            for k in range(i, n):
                sum = 0.0
                for j in range(i):
                    sum += (lower[i][j] * upper[j][k])
                upper[i][k] = matrix[i][k] - sum

            # Lower Triangular
            for k in range(i, n):
                if i == k:
                    lower[i][i] = 1.0
                else:
                    sum = 0.0
                    for j in range(i):
                        sum += (lower[k][j] * upper[j][i])
                    lower[k][i] = (matrix[k][i] - sum) / upper[i][i]

        return lower, upper

    lower, upper = _lu_decomposition(matrix)
    n = len(matrix)
    det = 1.0
    for i in range(n):
        det *= upper[i][i]
    return det


def identity(size: int):
    return [[1 if i == j else 0 for j in range(size)] for i in range(size)]


def durand_kerner(coefficients):
    n = len(coefficients) - 1  # degree of the polynomial
    roots = np.zeros(n, dtype=complex)  # array to store the roots

    # Initial guess for the roots
    for i in range(n):
        angle = 2.0 * np.pi * i / n
        roots[i] = np.cos(angle) + 1j * np.sin(angle)

    # Iteratively refine the roots
    max_iterations = 100
    tolerance = 1e-6
    for _ in range(max_iterations):
        delta = np.zeros(n, dtype=complex)
        for i in range(n):
            product = 1.0
            for j in range(n):
                if i != j:
                    product *= roots[i] - roots[j]
            delta[i] = -np.polyval(coefficients[:-1], roots[i]) / product
        roots += delta

        # Check convergence
        if np.all(np.abs(delta) < tolerance):
            break

    return roots


class DurandKernerRootsFinder:
    @staticmethod
    def _add_complex(c1, c2):
        return (c1[0] + c2[0], c1[1] + c2[1])

    @staticmethod
    def _multiply_complex(c1, c2):
        real = c1[0] * c2[0] - c1[1] * c2[1]
        imag = c1[0] * c2[1] + c1[1] * c2[0]
        return (real, imag)

    @staticmethod
    def _power_complex(c, n):
        result = (1, 0)
        for _ in range(n):
            result = DurandKernerRootsFinder._multiply_complex(result, c)
        return result

    @staticmethod
    def _evaluate_polynomial(a, c):
        n = len(a) - 1
        return sum([DurandKernerRootsFinder._multiply_complex(a[i], DurandKernerRootsFinder._power_complex(c, n - i))
                    for i in range(n + 1)], (0, 0))

    @staticmethod
    def calculate(a, eps=1e-6):
        n = len(a) - 1
        a = [(x[0] / a[0][0], x[1] / a[0][0]) for x in a]
        roots = [(1, i / (n + 1)) for i in range(n)]
        while True:
            new_roots = [
                DurandKernerRootsFinder._add_complex(
                    roots[i],
                    DurandKernerRootsFinder._multiply_complex((-1, 0), DurandKernerRootsFinder._multiply_complex(
                        DurandKernerRootsFinder._evaluate_polynomial(a, roots[i]),
                        (1 / abs(DurandKernerRootsFinder._evaluate_polynomial(a, roots[i])), 0))))
                for i in range(n)
            ]
            if max(abs(DurandKernerRootsFinder.add_complex(new_roots[i],
                                                           DurandKernerRootsFinder.multiply_complex((-1, 0), roots[i])))
                   for i in range(n)) < eps:
                return new_roots
            roots = new_roots


# class DeterminantTest(unittest.TestCase):
#     def test_random_size(self):
#         for size in range(1, 100):
#             matrix = np.random.random((size, size)).tolist()
#             np_det = np.linalg.det(matrix)
#             lu_det = lu_decomposition_determinant(matrix)
#
#             self.assertAlmostEqual(np_det, lu_det, places=1, msg=f'matrix size: {size}')
#

if __name__ == "__main__":
    matrix = np.random.random((20, 20)).tolist()

    roots = DurandKernerRootsFinder.calculate(matrix)
    print(roots)
