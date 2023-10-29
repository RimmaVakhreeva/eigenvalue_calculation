import cmath
import logging
import random
import unittest
from typing import List, Union, Tuple

from eigenval.utils import VectorType

__all__ = ['aberth_solve']


def _round_val(
        x: Union[float, complex],
        places: int = 6
) -> Union[Tuple[float, float], float]:
    """
    Round a number to a specified number of decimal places.

    If the number is complex, rounds both the real and imaginary parts separately.

    Parameters:
    - x (Union[float, complex]): The number (either real or complex) to be rounded.
    - places (int, optional): The number of decimal places to round to. Defaults to 6.

    Returns:
    - Union[Tuple[float, float], float]:
        If x is complex, returns a tuple of rounded real and imaginary parts.
        If x is a real number, returns the rounded value.
    """

    # Check if the input number is of complex type
    if type(x) is complex:
        # Round the real and imaginary parts separately and return as a tuple
        return round(x.real, places), round(x.imag, places)
    else:
        # If the number is not complex, simply round it
        return round(x, places)


def _aberth_method_trial(
        coefs: List[float],
        max_iter: int,
        tol: float,
        a: float,
        b: float
) -> List[Union[float, complex]]:
    """
    Approximates the roots of a polynomial using the Aberth method.

    Parameters:
    - coefs: List of polynomial coefficients, starting with the term of lowest degree.
    - max_iter: Maximum number of iterations for the Aberth method.
    - tol: Tolerance to determine convergence of the method.
    - a, b: Interval [a, b] for initial random root adjustments.

    Returns:
    - A list of roots, which can be either real or complex numbers.
    """

    # Helper function to evaluate a polynomial for a given x
    def _evaluate_polynomial(coefs, x):
        return sum(c * x ** i for i, c in enumerate(reversed(coefs)))

    # Helper function to evaluate the derivative of a polynomial for a given x
    def _evaluate_derivative(coefs, x):
        return sum(c * i * x ** (i - 1) for i, c in enumerate(reversed(coefs)) if i != 0)

    # Helper function to return real part of number if its imaginary part is negligibly small
    def _ignore_complex_if_small(num, threshold=1e-5):
        if abs(num.imag) < threshold:
            num = num.real
        return num

    # Determine the degree of the polynomial
    n = len(coefs) - 1

    # Initialize the roots with complex numbers distributed around the unit circle
    roots = [
        cmath.exp(2j * cmath.pi * i / n) * i * i
        for i in range(n)
    ]

    # Adjust each root with a random value between a and b
    for i in range(n):
        roots[i] += random.uniform(a, b)

    # Iteratively refine the roots using the Aberth method
    for _ in range(max_iter):
        new_roots = []
        for i in range(n):
            # Current approximation of the root for the polynomial
            root = roots[i]

            # Evaluate the polynomial using the current root approximation
            f = _evaluate_polynomial(coefs, root)

            # Evaluate the derivative of the polynomial using the current root approximation
            f_prime = _evaluate_derivative(coefs, root)

            # Calculate the differences between the current root and all other root approximations
            offsets = [root - other_root for other_root in roots]

            # Calculate the sum of reciprocal differences, excluding the one where offset is zero
            # (which corresponds to the root itself)
            sum_offsets = sum(1 / offset for offset in offsets if offset != 0)

            # Using Aberth's formula, compute the new approximation for the root.
            # The "+ 1e-10" in the denominator is to prevent division by zero or very small values,
            # which could cause instability.
            new_root = root - n * f / (n * f_prime - f * sum_offsets + 1e-10)

            # Store the new root approximation for the next iteration
            new_roots.append(new_root)

        # Check for convergence
        if max(abs(a - b) for a, b in zip(roots, new_roots)) < tol:
            return [_ignore_complex_if_small(num) for num in roots]

        # Update roots for next iteration
        roots = new_roots

    # Return refined roots after max_iter iterations
    return [_ignore_complex_if_small(num) for num in roots]


def aberth_solve(
        coefs: VectorType,
        trials: int = 100,
        iterations_per_trial: int = 100_000,
        tolerance: float = 1e-10,
        random_pert_lhs: float = -100.0,
        random_pert_rhs: float = 100.0,
) -> List[Union[float, complex]]:
    """
    Find roots of a polynomial using the Aberth method.

    The Aberth method is an iterative method used for finding all roots of a polynomial simultaneously.
    This function attempts to find these roots by performing multiple trials, each starting with a different
    set of initial approximations obtained through random perturbations. To ensure robustness against converging
    to the same root multiple times in different trials, roots are rounded and stored in a set (`seen`).

    Parameters:
    - coefs (VectorType): Coefficients of the polynomial in decreasing order of degree (
                          e.g., [a_n, a_{n-1}, ..., a_1, a_0] for a polynomial of degree n).
    - trials (int, optional): Number of trials to run the Aberth method. In each trial, a different initial estimate
                              is used. Default is 100.
    - iterations_per_trial (int, optional): Maximum number of iterations for each trial to attempt finding roots.
                                            Default is 10000.
    - tolerance (float, optional): Convergence criteria. If the difference between iterations is less than this value,
                                   it is assumed the method has converged. Default is 1e-10.
    - random_pert_lhs (float, optional): Left boundary for random perturbation of initial guesses. Default is -50.0.
    - random_pert_rhs (float, optional): Right boundary for random perturbation of initial guesses. Default is 50.0.

    Returns:
    - List[Union[float, complex]]: A list of unique roots of the polynomial. Depending on the polynomial,
                                   the roots can be real or complex numbers.

    Note:
    - The function depends on internal functions `_aberth_method_trial`
          (which performs the Aberth iteration for a given trial)
          and `_round_val` (which rounds a root for uniqueness checking)
    """
    if all(cmath.isclose(c, 0, rel_tol=1e-6, abs_tol=0.0) for c in coefs):
        return []

    # Initialize list to store unique roots and set to track seen roots
    unique_roots, seen = [], set()

    # Compute the degree of the polynomial from its coefficients
    n = len(coefs) - 1

    # Run Aberth method for multiple trials
    for i in range(trials):
        # Obtain roots for the current trial
        roots = _aberth_method_trial(
            coefs=coefs,
            max_iter=iterations_per_trial,
            tol=tolerance,
            a=random_pert_lhs,
            b=random_pert_rhs
        )

        # Iterate through each root found in the trial
        for root in roots:
            # Round the root for uniqueness checking
            rounded_root = _round_val(root, places=2)

            # If the rounded root has been seen before, skip
            if rounded_root in seen:
                continue

            # Add the rounded root to the seen set
            seen.add(rounded_root)

            # Append the actual root to the unique_roots list
            unique_roots.append(root)

        # Stop if the number of unique roots found matches the degree of the polynomial
        if len(unique_roots) == n:
            break

    # Return the list of unique roots
    if len(unique_roots) < n:
        logging.warning("Not all roots found")

    return unique_roots


class PolynomialRootsTest(unittest.TestCase):

    def test_linear_polynomial(self):
        coefs = [1, -5]  # x - 5
        roots = sorted(aberth_solve(coefs), key=lambda x: round(x.real, 4))
        self.assertAlmostEqual(roots[0], 5.0, places=4)

    def test_quadratic_real_roots(self):
        coefs = [1, -5, 6]  # x^2 - 5x + 6
        roots = sorted(aberth_solve(coefs), key=lambda x: round(x.real, 4))
        self.assertAlmostEqual(roots[0], 2.0, places=4)
        self.assertAlmostEqual(roots[1], 3.0, places=4)

    def test_quadratic_complex_roots(self):
        coefs = [1, 2, 5]  # x^2 + 2x + 5
        roots = sorted(aberth_solve(coefs), key=lambda x: (round(x.real, 4), round(x.imag, 4)))
        self.assertAlmostEqual(roots[0].real, -1.0, places=4)
        self.assertAlmostEqual(roots[0].imag, -2.0, places=4)
        self.assertAlmostEqual(roots[1].real, -1.0, places=4)
        self.assertAlmostEqual(roots[1].imag, 2.0, places=4)

    def test_all_zero_coefficients(self):
        coefs = [0, 0, 0]
        roots = aberth_solve(coefs)
        self.assertEquals(roots, [])

    def test_roots(self):
        import numpy as np

        for n in range(2, 20):
            coefs = (np.random.random(n) * 100 - 50).tolist()
            roots = sorted(aberth_solve(coefs), key=lambda x: (round(x.real, 4), round(x.imag, 4)))
            roots_np = sorted(np.roots(coefs).tolist(), key=lambda x: (round(x.real, 4), round(x.imag, 4)))

            self.assertEqual(len(roots), len(roots_np), msg=f"Failed coeff count: {n}")
            for root, root_np in zip(roots, roots_np):
                self.assertAlmostEqual(root, root_np, places=2, msg=f"Failed coeff count: {n}")


if __name__ == '__main__':
    unittest.main()
