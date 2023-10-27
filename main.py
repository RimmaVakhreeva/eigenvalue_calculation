def determinant_2x2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def determinant_3x3(matrix):
    return (matrix[0][0] * determinant_2x2([[matrix[1][1], matrix[1][2]],
                                            [matrix[2][1], matrix[2][2]]]) -
            matrix[0][1] * determinant_2x2([[matrix[1][0], matrix[1][2]],
                                            [matrix[2][0], matrix[2][2]]]) +
            matrix[0][2] * determinant_2x2([[matrix[1][0], matrix[1][1]],
                                            [matrix[2][0], matrix[2][1]]]))


def determinant_4x4(matrix):
    return (matrix[0][0] * determinant_3x3([row[1:] for row in matrix[1:]]) -
            matrix[0][1] * determinant_3x3([row[0:1] + row[2:] for row in matrix[1:]]) +
            matrix[0][2] * determinant_3x3([row[0:2] + row[3:] for row in matrix[1:]]) -
            matrix[0][3] * determinant_3x3([row[0:3] for row in matrix[1:]]))


def characteristic_polynomial(matrix, lmbda):
    n = len(matrix)
    identity = [[0] * n for _ in range(n)]
    for i in range(n):
        identity[i][i] = 1

    subtracted = [[matrix[i][j] - lmbda * identity[i][j] for j in range(n)] for i in range(n)]

    if n == 2:
        return determinant_2x2(subtracted)
    elif n == 3:
        return determinant_3x3(subtracted)
    elif n == 4:
        return determinant_4x4(subtracted)
    else:
        raise ValueError("Matrix size not supported.")


def find_eigenvalues(matrix):
    # Here, we're using the Newton-Raphson method to find the roots of the polynomial.
    # It's a simple numerical method but may not work for all matrices.
    # This implementation is very naive and for demonstration purposes.

    n = len(matrix)
    eigenvalues = []

    for _ in range(n):
        lmbda = 0  # Initial guess
        for _ in range(100):  # Iterate to refine the guess
            value = characteristic_polynomial(matrix, lmbda)
            derivative = (characteristic_polynomial(matrix, lmbda + 1e-6) - value) / 1e-6
            lmbda -= value / derivative  # Newton-Raphson update

        eigenvalues.append(lmbda)

    return eigenvalues


# Example
matrix = [[-4, -6], [3, 5]]
print(find_eigenvalues(matrix))



