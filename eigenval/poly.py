import unittest
from collections import defaultdict
from itertools import product
from random import randint
from typing import Dict

from eigenval.utils import MatrixType, VectorType, assert_matrix_is_square

__all__ = ['characteristic_polynomial']


def _minor(
        matrix: MatrixType,
        row: int,
        col: int
) -> MatrixType:
    """
    Compute the minor of a matrix after removing the specified row and column.

    Parameters:
    - matrix (MatrixType): A square matrix where each entry is either an int or a list of ints representing
                           polynomial coefficients.
    - row (int): The row index to be removed.
    - col (int): The column index to be removed.

    Returns:
    - MatrixType: A matrix with one fewer row and column, representing the minor of the original matrix.
    """
    return [
        [matrix[i][j] for j in range(len(matrix[i])) if j != col]
        for i in range(len(matrix)) if i != row
    ]


def _determinant_with_coeffs(matrix: MatrixType) -> Dict[int, float]:
    """
    Compute the determinant of a matrix, where each matrix element can be a polynomial.

    Given a matrix where each entry is either a scalar or a list of polynomial coefficients (with the index
    representing the degree of the polynomial term), this function computes the determinant and returns
    the coefficients of the resulting polynomial.

    Parameters:
    - matrix (MatrixType): A square matrix where each entry is either an int or a list of ints representing
                           polynomial coefficients.

    Returns:
    - Dict[int, float]: A dictionary mapping polynomial degrees to their respective coefficients. For example,
                        {0: 2, 1: -5, 2: 1} represents the polynomial 2 - 5*lambda + lambda^2.

    Note:
    This function uses the method of cofactor expansion, and it's recursive in nature. It's designed for matrices
    with polynomial entries, especially useful in computing characteristic polynomials.
    """
    # Get the number of rows (or columns) of the square matrix
    n = len(matrix)

    # Base case: If the matrix is 1x1, directly return the coefficients
    if n == 1:
        # Create a dictionary of coefficients from the single entry
        # in the matrix (which is a list of coefficients)
        return defaultdict(int, {i: coeff for i, coeff in enumerate(matrix[0][0])})

    # Dictionary to hold the accumulated coefficients for the polynomial of the determinant
    coefficients = defaultdict(int)

    # Loop through each column (or could be rows) for cofactor expansion
    for j in range(n):
        # Calculate the sign based on the column index.
        # Since we're expanding along the first row, only j affects the sign.
        sign = (-1) ** (j + 1)

        # Calculate the determinant (with polynomial coefficients)
        # of the minor matrix excluding the first row and j-th column
        minor_coeffs = _determinant_with_coeffs(_minor(matrix, 0, j))

        # Multiply the coefficients of the minor's determinant with the coefficients
        # of the matrix's element at (0, j)
        for (deg_m, coeff_m), (deg_cell, coeff_cell) in product(
                minor_coeffs.items(), enumerate(matrix[0][j])
        ):
            # Update the accumulated coefficient of the resulting polynomial
            # at the degree `deg_m + deg_cell`
            coefficients[deg_m + deg_cell] += sign * coeff_m * coeff_cell

    return coefficients


def characteristic_polynomial(matrix: MatrixType) -> VectorType:
    """
    Compute the characteristic polynomial of a square matrix.

    Given a square matrix `matrix`, this function computes its characteristic polynomial, which is defined as
    the determinant of `(matrix - lambda * I)`, where `lambda` is an indeterminate and `I` is the identity matrix.
    The result is returned as a list of coefficients in descending powers of `lambda`.

    Parameters:
    - matrix (MatrixType): A square matrix where each entry is either an int or a list of ints representing
                           polynomial coefficients. The matrix should be square, meaning its number of rows
                           should equal its number of columns.

    Returns:
    - VectorType: A list of coefficients representing the characteristic polynomial in descending powers of
                  `lambda`. For instance, the list [2, -5, 1] corresponds to the polynomial `2*lambda^2 - 5*lambda + 1`.

    Raises:
    - AssertionError: If the input matrix is not square.

    Note:
    The function utilizes cofactor expansion for computing the determinant with polynomial coefficients.
    """
    assert_matrix_is_square(matrix)

    n = len(matrix)
    if n == 0:
        return []

    # Create the matrix M - lambda * I and ensure all elements are lists
    augmented_matrix = [
        [[matrix[i][j]] if i != j else [matrix[i][j], -1] for j in range(n)]
        for i in range(n)
    ]

    coeffs_dict = _determinant_with_coeffs(augmented_matrix)
    coeffs_list = [-coeffs_dict[i] for i in range(max(coeffs_dict.keys()) + 1)]
    return list(reversed(coeffs_list))


class TestCharacteristicPolynomial(unittest.TestCase):

    def generate_random_matrix(self, size):
        """Generate a random matrix of the given size."""
        return [[randint(-10, 10) for _ in range(size)] for _ in range(size)]

    def assertListAlmostEqual(self, list1, list2, delta=1e-6):
        """Helper function to compare two lists with almost equal check."""
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, delta=delta)

    def test_on_random_matrices(self):
        import numpy as np

        """Test the characteristic polynomial on random matrices."""
        for n in range(1, 10):  # You can adjust the number of iterations
            matrix = self.generate_random_matrix(n)  # Adjust matrix size if needed
            our_poly = characteristic_polynomial(matrix)
            np_poly = list(np.poly(matrix))
            self.assertListAlmostEqual(our_poly, np_poly, delta=1e-4)

    def test_edge_case_empty_matrix(self):
        """Test on an empty matrix."""
        self.assertEqual(characteristic_polynomial([]), [])

    def test_edge_case_1x1_matrix(self):
        """Test on 1x1 matrix."""
        matrix = [[5]]
        self.assertEqual(characteristic_polynomial(matrix), [1, -5])


if __name__ == '__main__':
    unittest.main()
