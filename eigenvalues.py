import numpy as np

m = np.array([[5, 4, 2],
              [4, 5, 2],
              [2, 2, 2]])

w, v = np.linalg.eig(m)

print(w)


def determinant_2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def determinant_3(matrix):
    return (matrix[0][0] * determinant_2([[matrix[1][1], matrix[1][2],
                                         [matrix[2][1], matrix[2][2]]]]) -
            matrix[0][1] * determinant_2(([[matrix[1][0], matrix[1][2],
                                          [matrix[2][0], matrix[2][2]]]]) +
            matrix[0][2] * determinant_2([[matrix[1][0], matrix[1][1]],
                                         [matrix[2][0], matrix[2][1]]]))
            )


def determinant_4(matrix):
    return (matrix[0][0] * determinant_3([row[1:] for row in matrix[1:]]) -
            matrix[0][1] * determinant_3([row[0:1] + row[2:] for row in matrix[1:]]) +
            matrix[0][2] * determinant_3([row[0:2] + row[3:] for row in matrix[1:]]) -
            matrix[0][3] * determinant_3([row[0:3] for row in matrix[1:]])
            )


def characteristic_equation(matrix, eigenvalue_lambda):
    n = len(matrix)
    identity = [[0] * n for _ in range(n)]
    for i in range(n):
        identity[i][i] = 1

    sub_eigenvalue_lambda = [[matrix[i][j] - eigenvalue_lambda * identity[i][j] for j in range(n)] for i in range(n)]

    if n == 2:
        return determinant_2(sub_eigenvalue_lambda)
    if n == 3:
        return determinant_3(sub_eigenvalue_lambda)
    if n == 4:
        return determinant_4(sub_eigenvalue_lambda)
    else:
        raise ValueError('Matrix size is not supported')


def find_eigenvalues(matrix):
    n = len(matrix)
    epsilon = 0.000001

    for _ in range(n):
        eigenvalue_0 = -5
        eigenvalue_1 = 5
        idx = 0
        while abs(eigenvalue_1 - eigenvalue_0) > epsilon:
            tmp = eigenvalue_1
            eigenvalue_1 = eigenvalue_1 - characteristic_equation(matrix, eigenvalue_1) * ((eigenvalue_1 - eigenvalue_0) / (characteristic_equation(matrix, eigenvalue_1) - characteristic_equation(matrix, eigenvalue_0)))
            eigenvalue_0 = tmp
            idx += 1

        return eigenvalue_0, eigenvalue_1


matrix = [[5, 4, 2],
          [4, 5, 2],
          [2, 2, 2]]
print(find_eigenvalues(matrix))

