import cmath
import math

_math_functions = [
    math.cos,
    math.sin,
]


def linear_function(x):
    return x


def quadratic_function(x):
    return x ** 2


def poly_function(x):
    return x ** 3 + 2 * x ** 2 - 5 * x + 1


def custom_func_1(z):
    return z ** 2 + 3 * z + 2


def custom_func_2(z):
    return math.exp(z) + z ** 3 - 2


def custom_func_3(z):
    return math.sin(z) + math.cos(z)


def custom_func_4(z):
    return (z + 1) / (z - 2)


def custom_func_5(z):
    return math.exp(z) * math.sin(z)


def custom_func_6(z):
    return abs(z) * cmath.phase(z)


all_functions = [
    *_math_functions,
    linear_function,
    quadratic_function,
    poly_function,
    custom_func_1,
    custom_func_2,
    custom_func_3,
    custom_func_5,
    custom_func_6,
]
