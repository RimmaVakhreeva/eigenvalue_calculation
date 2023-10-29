import unittest

from eigenval import roots, poly
from eigenval.utils import MatrixType, VectorType

__all__ = ['eigvals']


def eigvals(matrix: MatrixType) -> VectorType:
    """
    Compute the eigenvalues of a given matrix using its characteristic polynomial.

    Parameters:
    - matrix (MatrixType): Input matrix for which the eigenvalues are to be calculated.

    Returns:
    - VectorType: Vector containing the eigenvalues of the matrix.
    """

    # Compute the characteristic polynomial of the matrix
    char_poly = poly.characteristic_polynomial(matrix)

    # Find the roots of the characteristic polynomial to get the eigenvalues
    return roots.aberth_solve(char_poly)


class TestEigvals(unittest.TestCase):

    def test_eigvals(self):
        import numpy as np

        for n in range(1, 10):
            matrix = np.random.rand(n, n)

            custom_eigvals = eigvals(matrix.tolist())
            numpy_eigvals = np.linalg.eigvals(matrix)

            custom_eigvals = sorted(custom_eigvals, key=lambda x: (round(x.real, 4), round(x.imag, 4)))
            numpy_eigvals = sorted(numpy_eigvals.tolist(), key=lambda x: (round(x.real, 4), round(x.imag, 4)))

            self.assertEqual(len(custom_eigvals), len(numpy_eigvals), msg=f"Failed val count: {n}")
            for custom_val, np_val in zip(custom_eigvals, numpy_eigvals):
                self.assertAlmostEqual(custom_val, np_val, places=2, msg=f"Failed matrix size: {n}")


if __name__ == "__main__":
    unittest.main()
