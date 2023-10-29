from typing import List

MatrixType = List[List[float]]
VectorType = List[float]


def assert_matrix_is_square(matrix: MatrixType):
    n = len(matrix)
    for row in matrix:
        assert len(row) == n, f"Matrix {matrix} is not square"
