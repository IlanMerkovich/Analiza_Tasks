def lagrange_interpolation(x_data, y_data, x):
    """
    Lagrange Interpolation

    Parameters:
    x_data (list): List of x-values for data points.
    y_data (list): List of y-values for data points.
    x (float): The x-value where you want to evaluate the interpolated polynomial.

    Returns:
    float: The interpolated y-value at the given x.
    """
    n = len(x_data)
    result = 0.0

    for i in range(n):
        term = y_data[i]
        for j in range(n):
            if i != j:
                term *= (x - x_data[j]) / (x_data[i] - x_data[j])
        result += term

    return result


def neville(x_data, y_data, x_interpolate):
    """
    Neville's Interpolation

    Parameters:
    x_data (list): List of x-values for data points.
    y_data (list): List of y-values for data points.
    x_interpolate (float): The x-value where you want to evaluate the interpolated polynomial.

    Returns:
    float: The interpolated y-value at the given x_interpolate.
    """
    n = len(x_data)
    tableau = [[0.0] * n for _ in range(n)]
    for i in range(n):
        tableau[i][0] = y_data[i]
    for j in range(1, n):
        for i in range(n - j):
            tableau[i][j] = ((x_interpolate - x_data[i + j]) * tableau[i][j - 1] -
                             (x_interpolate - x_data[i]) * tableau[i + 1][j - 1]) / \
                            (x_data[i] - x_data[i + j])

    return tableau[0][n - 1]


if __name__ == "__main__":
    x_data = [1, 2, 3, 4]
    y_data = [1, 4, 9, 16]
    x_point = 2.5

    lagrange_result = lagrange_interpolation(x_data, y_data, x_point)
    print(f"Lagrange Interpolation Result at x = {x_point}: {lagrange_result}")

    neville_result = neville(x_data, y_data, x_point)
    print(f"Neville Interpolation Result at x = {x_point}: {neville_result}")
