import random
import unittest
from typing import Callable

from integrate.funcs import all_functions

__all__ = ['trapezoid']


def trapezoid(func: Callable[[float], float], a: float, b: float, n: int) -> float:
    """
    Approximate the integral of a function using the trapezoidal rule.

    Parameters:
    - func: A callable function which takes a float and returns a float.
            Represents the function to be integrated.
    - a: The starting point of the interval over which to integrate.
    - b: The ending point of the interval over which to integrate.
    - n: The number of trapezoids or subdivisions.

    Returns:
    - The approximated integral value as a float.

    Note:
    The function assumes that `a` is less than or equal to `b` and `n` is positive.
    """

    # Check prerequisites
    assert callable(func), "func should be a callable function"
    assert n > 0, "Number of trapezoids (n) should be positive"
    assert a <= b, "Start of interval (a) should be less than or equal to end of interval (b)"

    # Calculate the width of each trapezoid
    h = float((b - a) / n)

    # Compute the sum of the first and last y-values
    sum = 0.5 * (func(a) + func(b))

    # Compute the sum of the y-values for the interior points
    for i in range(1, n):
        sum += func(a + i * h)

    # Multiply by the width of the trapezoids to get the final integral approximation
    return sum * h


class TestTrapezoid(unittest.TestCase):

    def test_accuracy_against_scipy(self):
        from scipy.integrate import trapezoid as scipy_trapezoid
        import numpy as np

        for func in all_functions:
            a = random.uniform(-1, 1)
            b = random.uniform(a + 1e-2, a + 1)
            n = random.randint(10, 1000)

            with self.subTest(func=func):
                x = np.linspace(a, b, num=n).tolist()
                y = np.array([func(xi) for xi in x], dtype=np.float32)

                try:
                    custom_val = trapezoid(func, a, b, n)
                    scipy_val = scipy_trapezoid(y, x)
                    self.assertAlmostEqual(
                        custom_val,
                        scipy_val,
                        delta=1e-3,
                        msg=f"func: {func}"
                    )
                except ValueError as e:
                    print(f"func: {func}")
                    raise e

    def test_negative_n(self):
        func = all_functions[0]  # Just taking the first function for simplicity
        with self.assertRaises(AssertionError):
            trapezoid(func, 0, 2, -100)

    def test_a_greater_than_b(self):
        func = all_functions[0]  # Just taking the first function for simplicity
        with self.assertRaises(AssertionError):
            trapezoid(func, 2, 0, 100)

    def test_n_as_zero(self):
        func = all_functions[0]  # Just taking the first function for simplicity
        with self.assertRaises(AssertionError):
            trapezoid(func, 0, 2, 0)


if __name__ == "__main__":
    unittest.main()
