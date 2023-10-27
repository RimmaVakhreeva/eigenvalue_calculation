import cmath
import numpy as np
import unittest


def aberth_method(coefs, max_iter=1000, tol=1e-10):
    def evaluate_polynomial(coefs, x):
        return sum(c * x ** i for i, c in enumerate(reversed(coefs)))

    def evaluate_derivative(coefs, x):
        return sum(c * i * x ** (i - 1) for i, c in enumerate(reversed(coefs)) if i != 0)

    def _ignore_complex_if_small(num, threshold=1e-5):
        if abs(num.imag) < threshold:
            num = num.real
        return num

    n = len(coefs) - 1
    roots = [
        cmath.exp(2j * cmath.pi * i / n) * i*i
        for i in range(n)
    ]
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
            return [_ignore_complex_if_small(num) for num in roots]
        roots = new_roots

    return [_ignore_complex_if_small(num) for num in roots]


class PolynomialRootsTest(unittest.TestCase):
    def test_roots(self):
        for i in range(3, 100):
            coefs = np.random.random(i).tolist()
            roots = sorted(aberth_method(coefs), key=lambda x: (round(x.real, 4), round(x.imag, 4)))
            roots_np = sorted(np.roots(coefs).tolist(), key=lambda x: (round(x.real, 4), round(x.imag, 4)))

            for root, root_np in zip(roots, roots_np):
                self.assertAlmostEqual(root, root_np, places=4, msg=f"{i} iter")
