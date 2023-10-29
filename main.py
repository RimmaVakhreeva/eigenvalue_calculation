import math

from eigenval import eigvals
from integrate import trapezoid


def task_1():
    def _pretty_print(m):
        message = ''
        for row in m:
            message += ' '.join(map(str, row))
            message += '\n'
        return message.strip()

    matrix = [
        [0.12, 0.45, 0.73, 0.56, 0.89],
        [0.65, 0.23, 0.91, 0.38, 0.74],
        [0.47, 0.85, 0.62, 0.29, 0.18],
        [0.39, 0.14, 0.97, 0.53, 0.81],
        [0.69, 0.77, 0.36, 0.11, 0.57]
    ]
    vals = eigvals(matrix)
    print(f"Task 1. Given the matrix:\n{_pretty_print(matrix)}")
    print(f"Eigenvalues are: {[round(v, 2) for v in vals]}\n")


def task_2():
    def _func(x):
        return math.exp(x) + x ** 3 - 2

    integral_value = trapezoid(
        func=_func,
        a=-1,
        b=1,
        n=100
    )
    print(f"Task 2. Given the function f(x) = e^x + x^3 - 2")
    print(f"Integral value from -1 to 1 is: {integral_value:.2f}")


if __name__ == '__main__':
    task_1()
    task_2()
