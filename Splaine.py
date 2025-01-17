import numpy as np

def spline_cubic(x_list, y_list, x, f_tag_0, f_tag_n):
    """
    Cubic Spline Interpolation with Natural and Full Spline.

    Parameters:
    x_list (list): List of x-values for data points.
    y_list (list): List of y-values for data points.
    x (float): The x-value to predict.
    f_tag_0 (float): Derivative at the first point (for full spline).
    f_tag_n (float): Derivative at the last point (for full spline).

    Returns:
    tuple: (natural_spline_result, full_spline_result)
    """
    n = len(x_list) - 1
    h_list = [x_list[i+1] - x_list[i] for i in range(n)]

    # Step 1: Build lambda, mu, and d for the natural spline
    lam_list = [h_list[i] / (h_list[i] + h_list[i-1]) if i > 0 else 0 for i in range(n)]
    mu_list = [1 - lam_list[i] if i > 0 else 0 for i in range(n)]
    d_list = [0] * (n + 1)
    for i in range(1, n):
        d_list[i] = 6 * ((y_list[i+1] - y_list[i]) / h_list[i] -
                         (y_list[i] - y_list[i-1]) / h_list[i-1]) / (h_list[i] + h_list[i-1])

    # Step 2: Solve the system of equations for natural spline
    matrix = np.zeros((n+1, n+1))
    for i in range(1, n):
        matrix[i][i-1] = mu_list[i]
        matrix[i][i] = 2
        matrix[i][i+1] = lam_list[i]
    matrix[0][0] = 1
    matrix[n][n] = 1
    m_natural = np.linalg.solve(matrix, d_list)

    # Step 3: Calculate the result for the natural spline
    for i in range(n):
        if x_list[i] <= x <= x_list[i+1]:
            dx = x - x_list[i]
            result_natural = (y_list[i] +
                              (m_natural[i] / 2) * dx**2 +
                              ((m_natural[i+1] - m_natural[i]) / (6 * h_list[i])) * dx**3 +
                              ((y_list[i+1] - y_list[i]) / h_list[i] - (h_list[i] / 6) *
                               (m_natural[i+1] + 2 * m_natural[i])) * dx)
            break

    # Step 4: Build the matrix and vector for the full spline
    d_list[0] = 6 * ((y_list[1] - y_list[0]) / h_list[0] - f_tag_0) / h_list[0]
    d_list[n] = 6 * (f_tag_n - (y_list[n] - y_list[n-1]) / h_list[n-1]) / h_list[n-1]

    matrix[0][0] = 2
    matrix[0][1] = 1
    matrix[n][n-1] = 1
    matrix[n][n] = 2
    m_full = np.linalg.solve(matrix, d_list)

    # Step 5: Calculate the result for the full spline
    for i in range(n):
        if x_list[i] <= x <= x_list[i+1]:
            dx = x - x_list[i]
            result_full = (y_list[i] +
                           (m_full[i] / 2) * dx**2 +
                           ((m_full[i+1] - m_full[i]) / (6 * h_list[i])) * dx**3 +
                           ((y_list[i+1] - y_list[i]) / h_list[i] - (h_list[i] / 6) *
                            (m_full[i+1] + 2 * m_full[i])) * dx)
            break

    return result_natural, result_full


if __name__ == "__main__":
    x_list = [0, 1, 2, 3]
    y_list = [1, 2, 0, 2]
    x = 1.5
    f_tag_0 = 1.0
    f_tag_n = -1.0
    natural_result, full_result = spline_cubic(x_list, y_list, x, f_tag_0, f_tag_n)

    print(f"Natural Spline Result at x = {x}: {natural_result:.4f}")
    print(f"Full Spline Result at x = {x}: {full_result:.4f}")
