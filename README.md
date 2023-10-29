# Assesment 2

This project contains code files related to eigenvalues and integration tasks. It includes functions to calculate eigenvalues of a matrix, compute the characteristic polynomial of a matrix, find roots of a polynomial, and approximate the integral of a given function using the trapezoidal rule.

## Project Files

The project folder structure is as follows:

```
- main.py
-- eigenval/eigenvalues.py
-- eigenval/poly.py
-- eigenval/roots.py
-- integrate/funcs.py
-- integrate/impl.py
```

The project files and their descriptions are listed below:

- `main.py`: This code file contains two tasks: `task_1` and `task_2`. In `task_1`, a matrix is defined, and its eigenvalues are calculated using the `eigvals()` function from the `eigenval` module. The matrix and calculated eigenvalues are then printed. In `task_2`, a function `f(x) = e^x + x^3 - 2` is defined. The integral value of this function from -1 to 1 is calculated using the `trapezoid()` function from the `integrate` module. The calculated integral value is then printed.

### eigenval

The `eigenval` module contains code files related to the calculation of eigenvalues of a matrix.

- `eigenvalues.py`: This code file contains a function called `eigvals` that calculates the eigenvalues of a given matrix using its characteristic polynomial. The function takes a matrix as input and returns a vector containing the eigenvalues. It computes the characteristic polynomial of the matrix using the `characteristic_polynomial` function from the `poly` module and then uses the `aberth_solve` function from the `roots` module to find the roots of the characteristic polynomial, which are the eigenvalues of the matrix. The code also includes a unit test class called `TestEigvals` that tests the `eigvals` function. It uses the `numpy` library to generate random matrices of different sizes and compares the eigenvalues calculated by the `eigvals` function with the eigenvalues calculated by the `numpy.linalg.eigvals` function. The test asserts that the number of eigenvalues is the same and checks if each eigenvalue is approximately equal within a certain number of decimal places.

- `poly.py`: The `poly` module contains a function called `characteristic_polynomial` that computes the characteristic polynomial of a square matrix. The characteristic polynomial is defined as the determinant of `(matrix - lambda * I)`, where `matrix` is the input matrix, `lambda` is an indeterminate, and `I` is the identity matrix. The function returns a list of coefficients representing the characteristic polynomial in descending powers of `lambda`. The code also includes a helper function called `_determinant_with_coeffs` that calculates the determinant of a matrix with polynomial entries using cofactor expansion. This function is used internally by `characteristic_polynomial`. The code includes a unit test class called `TestCharacteristicPolynomial` that tests the `characteristic_polynomial` function on random matrices. It also includes additional edge case tests for an empty matrix and a 1x1 matrix.

- `roots.py`: The `roots` module provides an implementation of the Aberth method for finding the roots of a polynomial. The main function is `aberth_solve`, which takes the coefficients of the polynomial as input and returns a list of unique roots as output. The Aberth method is an iterative method that uses initial approximations obtained through random perturbations to find all roots of a polynomial simultaneously. The method is run for multiple trials with different initial estimates to ensure robustness against converging to the same root multiple times. The code also includes helper functions such as `_round_val` to round numbers to a specified number of decimal places and `_aberth_method_trial` to perform the Aberth iteration for a given trial. The code includes a test case class called `PolynomialRootsTest` that verifies the correctness of the implementation using various test cases. The test cases cover linear polynomials, quadratic polynomials with real and complex roots, polynomials with all zero coefficients, and randomly generated polynomials with varying degrees.

### integrate

The `integrate` module contains code files related to the approximation of integrals using the trapezoidal rule.

- `funcs.py`: This code file defines several mathematical functions, including linear, quadratic, and polynomial functions. It also includes custom functions that involve complex numbers and trigonometric functions. All the functions are stored in a list called `all_functions`.

- `impl.py`: This code file contains a function called `trapezoid` that approximates the integral of a given function using the trapezoidal rule. It takes as input a callable function, the starting and ending points of the interval, and the number of trapezoids or subdivisions. The function first checks that the input meets certain prerequisites, such as ensuring that the function is callable, the number of trapezoids is positive, and the start of the interval is less than or equal to the end of the interval. Next, it calculates the width of each trapezoid and initializes a variable to store the sum of the y-values. It then computes the sum of the first and last y-values and iterates over the interior points to compute their sum as well. Finally, it multiplies the sum by the width of the trapezoids to obtain the final integral approximation, which is returned as the output. The code also includes a test class called `TestTrapezoid`, which contains several test methods. The `test_accuracy_against_scipy` method compares the result of the `trapezoid` function with the result of the `scipy.trapezoid` function for various functions and inputs. The `test_negative_n`, `test_a_greater_than_b`, and `test_n_as_zero` methods test the behavior of the `trapezoid` function for invalid inputs.

## How to Run

To run the code, follow the instructions below:

1. Execute the `main.py` script.

## Conclusion

This project provides functions for calculating eigenvalues of a matrix, computing the characteristic polynomial of a matrix, finding roots of a polynomial, and approximating integrals using the trapezoidal rule. Each code file includes relevant documentation and unit tests to ensure the correctness of the implementations.