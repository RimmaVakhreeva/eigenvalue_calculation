def f(x):
    # Define your function here. For example:
    return x ** 2


def trapezoidal_integral(func, a, b, n):
    """
    Estimates the integral of a function using the trapezoidal rule.

    Parameters:
    - func: the function to be integrated.
    - a: the start of the integration interval.
    - b: the end of the integration interval.
    - n: number of trapezoids.

    Returns:
    - Estimated integral value.
    """
    h = (b - a) / n
    sum = 0.5 * (func(a) + func(b))

    for i in range(1, n):
        sum += func(a + i * h)

    return sum * h


if __name__ == "__main__":
    # Define the interval [a, b] and number of trapezoids
    a = 0
    b = 3
    n = 1000

    result = trapezoidal_integral(f, a, b, n)
    print(f"Estimated integral value: {result:.4f}")

