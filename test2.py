import cmath

import numpy as np


def determinant(matrix):
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


def characteristic_polynomial(matrix):
    n = len(matrix)
    if n != len(matrix[0]):
        raise ValueError("Input matrix must be square")

    def cofactor(matrix, i, j):
        sub_matrix = [[matrix[row][col] for col in range(n) if col != j] for row in range(n) if row != i]
        return ((-1) ** (i + j)) * determinant(sub_matrix)

    def adjugate(matrix):
        return [[cofactor(matrix, i, j) for j in range(n)] for i in range(n)]

    char_poly = [1]
    adj_matrix = adjugate(matrix)
    for k in range(1, n + 1):
        char_poly.append((-1) ** k * sum(adj_matrix[k - 1]))

    return char_poly


def aberth_method(coefs, max_iter=100, tol=1e-10):
    def evaluate_polynomial(coefs, x):
        return sum(c * x ** i for i, c in enumerate(reversed(coefs)))

    def evaluate_derivative(coefs, x):
        return sum(c * i * x ** (i - 1) for i, c in enumerate(reversed(coefs)) if i != 0)

    n = len(coefs) - 1
    roots = [cmath.exp(2j * cmath.pi * i / n) for i in range(n)]
    for _ in range(max_iter):
        new_roots = []
        for i in range(n):
            root = roots[i]
            f = evaluate_polynomial(coefs, root)
            f_prime = evaluate_derivative(coefs, root)
            offsets = [root - other_root for other_root in roots]
            sum_offsets = sum(1 / offset for offset in offsets if offset != 0)
            new_root = root - n * f / (n * f_prime - f * sum_offsets + 1e-10)
            new_roots.append(new_root)
        if max(abs(a - b) for a, b in zip(roots, new_roots)) < tol:
            return new_roots
        roots = new_roots
    return roots


def eigenvalues(matrix):
    char_poly = characteristic_polynomial(matrix)
    roots = aberth_method(char_poly)
    roots_np = np.roots(char_poly)
    return roots


if __name__ == "__main__":
    matrix = np.random.random((2, 2)).tolist()
    ev = eigenvalues(matrix)
    ev_np, _ = np.linalg.eig(matrix)
    ev_np = ev_np.tolist()
    print(ev)
