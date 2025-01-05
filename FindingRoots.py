import sympy as sp
def find_intervals(func, start, end, step=0.1):
    """
    Finds intervals [x1, x2] where f(x1)*f(x2) <= 0, indicating a root exists between x1 and x2.
   """
    intervals = []
    x1 = start
    f_x1 = func(x1)

    while x1 < end:
        x2 = x1 + step
        f_x2 = func(x2)

        if f_x1 * f_x2 <= 0:
            intervals.append((x1, x2))

        x1 = x2
        f_x1 = f_x2

    print(intervals)
    return intervals


def bisection_method(func, start_point, end_point, tolerance):
    """
    Implements the Bisection Method to find a root of a function.
    """
    iterations = 0
    while (end_point - start_point) / 2 > tolerance:
        iterations += 1
        midpoint = (start_point + end_point) / 2
        if func(midpoint) == 0:
            return midpoint, iterations
        elif func(start_point) * func(midpoint) < 0:
            end_point = midpoint
        else:
            start_point = midpoint
    root = (start_point + end_point) / 2
    return root, iterations

def secant_method(func, start_point, end_point, tolerance, max_iterations=1000):
    """
    Implements the Secant Method to find a root of a function.
    """
    x0, x1 = start_point, end_point
    iterations = 0

    for _ in range(max_iterations):
        iterations += 1
        f_x0 = func(x0)
        f_x1 = func(x1)

        if f_x1 - f_x0 == 0:
            print("Error: Division by zero.")
            return None, iterations

        x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)

        if abs(x2 - x1) < tolerance:
            return x2, iterations

        x0, x1 = x1, x2

    print("method did not converge.")
    return None, iterations


def newton_raphson_method(func, derivative, initial_guess, tolerance=0.0001, max_iterations=1000):
    """
    Implements the Newton-Raphson Method to find a root of a function.
    """
    x_current = initial_guess
    iterations = 0

    for _ in range(max_iterations):
        iterations += 1
        f_value = func(x_current)
        derivative_value = derivative(x_current)

        if derivative_value == 0:
            print("Error: Derivative is zero.")
            return None, iterations

        x_next = x_current - f_value / derivative_value

        if abs(x_next - x_current) < tolerance:
            return x_next, iterations

        x_current = x_next

    print("method did not converge.")
    return None, iterations


if __name__ == "__main__":
    def sample_function(x):
        return x**3-4*x**2+3

    def sample_derivative(x):
        return 3*x**2-8*x

    start = -10
    end = 10
    step = 0.1
    tolerance = 0.0000001

    print("Finding intervals where roots exist:")
    intervals = find_intervals(sample_function, start, end, step)
    print("Intervals:", intervals)

    print("\nFinding roots using Bisection Method:")
    for interval in intervals:
        root, iterations = bisection_method(sample_function, interval[0], interval[1], tolerance)
        if root is not None:
            print(f"Root found using Bisection: {root} (in {iterations} iterations)")

    print("\nFinding roots using Secant Method:")
    for interval in intervals:
        root, iterations = secant_method(sample_function, interval[0], interval[1], tolerance)
        if root is not None:
            print(f"Root found using Secant: {root} (in {iterations} iterations)")

    print("\nFinding roots using Newton-Raphson Method:")
    for interval in intervals:
        initial_guess = (interval[0] + interval[1]) / 2
        root, iterations = newton_raphson_method(sample_function, sample_derivative, initial_guess, tolerance)
        if root is not None:
            print(f"Root found using Newton-Raphson: {root} (in {iterations} iterations)")
