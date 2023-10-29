The code file contains a function `characteristic_polynomial` that computes the characteristic polynomial of a square matrix. The characteristic polynomial is defined as the determinant of `(matrix - lambda * I)`, where `matrix` is the input matrix, `lambda` is an indeterminate, and `I` is the identity matrix. The function returns a list of coefficients representing the characteristic polynomial in descending powers of `lambda`.

The code also includes a helper function `_determinant_with_coeffs` that calculates the determinant of a matrix with polynomial entries using cofactor expansion. This function is used internally by `characteristic_polynomial`.

The code includes a unit test class `TestCharacteristicPolynomial` that tests the `characteristic_polynomial` function on random matrices. It also includes additional edge case tests for an empty matrix and a 1x1 matrix.

To run the code and execute the tests, the script can be executed directly.