import math
import sympy as sp


def max_iterations(a, b, err):
    s = int(math.floor(-math.log2(err / (b - a)) - 1))
    return s


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

    return intervals


def bisection_method(f, a, b, tol=1e-6):
    """
    Perform the Bisection Method to find a root of f in [a, b].
    """
    if math.copysign(1, f(a)) == math.copysign(1, f(b)):
        raise Exception("The scalars a and b do not bound a root")

    c, k = 0, 0
    max_iter = max_iterations(a, b, tol)

    while abs(b - a) > tol and k < max_iter:
        c = a + (b - a) / 2

        if f(c) == 0:
            return c

        if f(c) * f(a) < 0:
            b = c
        else:
            a = c
        k += 1

    return c,k

def secant_method(f, x0, x1, TOL, N=100):
    """
    Perform the Secant Method to find a root of f.
    """
    for i in range(N):
        if f(x1) - f(x0) == 0:
            print("Error: Division by zero in Secant Method.")
            return None, i

        p = x0 - f(x0) * ((x1 - x0) / (f(x1) - f(x0)))

        if abs(p - x1) < TOL:
            return p, i + 1

        x0 = x1
        x1 = p

    print("Secant Method did not converge within the maximum iterations.")
    return None, N


def newton_raphson_method(f, df, p0, TOL, N=100):
    """
    Perform the Newton-Raphson Method to find a root of f starting at p0.
    """
    for i in range(N):
        if df(p0) == 0:
            print("Error: Derivative is zero at p0 in Newton-Raphson Method.")
            return None, i

        p = p0 - f(p0) / df(p0)

        if abs(p - p0) < TOL:
            return p, i + 1

        p0 = p

    print("Newton-Raphson Method did not converge within the maximum iterations.")
    return None, N


if __name__ == "__main__":
    def sample_function(x):
        return x**2 - 4 * math.sin(x)


    def sample_derivative(x):
        return 3*x**2 - 4*math.cos(x)

    start = -100
    end = 100
    step = 0.1
    tolerance = 1e-16

    print("Finding intervals where roots exist:")
    intervals = find_intervals(sample_function, start, end, step)
    if intervals:
        print("Intervals:", intervals)
        print("\nFinding roots using Bisection Method:")
        for interval in intervals:
            try:
                root, iterations = bisection_method(sample_function, interval[0], interval[1], tolerance)
                if root is not None:
                    print(f"Root found using Bisection: {root:.10f} (in {iterations} iterations)")
            except Exception as e:
                print(f"Bisection Method Error: {e}")

        print("\nFinding roots using Secant Method:")
        for interval in intervals:
            root, iterations = secant_method(sample_function, interval[0], interval[1], tolerance)
            if root is not None:
                print(f"Root found using Secant: {root:.10f} (in {iterations} iterations)")

        print("\nFinding roots using Newton-Raphson Method:")
        for interval in intervals:
            initial_guess = (interval[0] + interval[1]) / 2
            root, iterations = newton_raphson_method(sample_function, sample_derivative, initial_guess, tolerance)
            if root is not None:
                print(f"Root found using Newton-Raphson: {root:.10f} (in {iterations} iterations)")
    else:
        print("No roots found!")
