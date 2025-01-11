
def secant_method(f, x0, x1, TOL=1e-10, N=50):
    for i in range(N):
        if f(x1) - f(x0) == 0:
            return None, i
        p = x0 - f(x0) * ((x1 - x0) / (f(x1) - f(x0)))
        if abs(p - x1) < TOL:
            return p, i + 1
        x0 = x1
        x1 = p
    return None, N


# Newton-Raphson Method
def newton_raphson(f, df, p0, TOL=1e-10, N=50):
    for i in range(N):
        if df(p0) == 0:
            return None, i
        p = p0 - f(p0) / df(p0)
        if abs(p - p0) < TOL:
            return p, i + 1
        p0 = p
    return None, N


# Bisection Method
def bisection_method(f, a, b, tol=1e-10, max_steps=50):
    if f(a) * f(b) >= 0:
        return None, 0
    c = 0
    for k in range(max_steps):
        c = (a + b) / 2  # Midpoint
        if f(c) == 0 or abs(b - a) < tol:
            return c, k + 1
        if f(c) * f(a) < 0:
            b = c  # Root is in [a, c]
        else:
            a = c  # Root is in [c, b]
    return c, max_steps


# Find intervals with potential roots
def find_potential_root_intervals(f, df, start, end, step=0.1):
    intervals = []
    x = start
    while x < end:
        x_next = x + step
        if f(x) * f(x_next) < 0:
            intervals.append((x, x_next))
        elif f(x) == 0 or df(x) * df(x_next) < 0:
            intervals.append((x, x_next))
        x = x_next
    return intervals


# Main Program
def main():
    print("=" * 50)
    print("Running Root-Finding Program")
    print("=" * 50)

    # Define the function and its derivative
    def f(x):
        return x ** 4 + x ** 3 - 3 * x ** 2  # Example function

    def df(x):
        return 4 * x ** 3 + 3 * x ** 2 - 6 * x  # Derivative of the function

    # Define the range and step size
    start, end, step = -10, 10, 0.1

    # Find intervals and roots
    intervals = find_potential_root_intervals(f, df, start, end, step)

    print("\nIntervals with potential simple or flat roots in f(x):")
    print("-" * 50)
    for interval in intervals:
        print(f"Interval: {interval[0]:.6f} to {interval[1]:.6f}")
    print("-" * 50)

    print("\nChoose a method to solve for the roots.\n")

    while True:
        print("=" * 50)
        print("Choose a method:")
        print("1. Secant Method")
        print("2. Newton-Raphson Method")
        print("3. Bisection Method")
        print("0. Exit")
        print("=" * 50)

        choice = input("Enter the number of your choice: ")

        if choice == '0':
            print("\nExiting the program. Goodbye!")
            break

        # Run the chosen method for each interval
        for interval in intervals:
            a, b = interval
            print(f"\nProcessing interval [{a:.6f}, {b:.6f}]:")

            if choice == '1':
                # Secant Method
                secant_root, secant_iterations = secant_method(f, a, b)
                if secant_root is not None:
                    print(f"Secant method root: {secant_root:.6f}")
                    print(f"Number of iterations: {secant_iterations}")
                else:
                    print("Secant method failed.")

            elif choice == '2':
                # Newton-Raphson Method (using midpoint of interval as starting point)
                midpoint = (a + b) / 2
                newton_root, newton_iterations = newton_raphson(f, df, midpoint)
                if newton_root is not None:
                    print(f"Newton-Raphson method root: {newton_root:.6f}")
                    print(f"Number of iterations: {newton_iterations}")
                else:
                    print("Newton-Raphson method failed.")

            elif choice == '3':
                # Bisection Method
                try:
                    bisection_root, bisection_iterations = bisection_method(f, a, b)
                    if bisection_root is not None:
                        print(f"Bisection method root: {bisection_root:.6f}")
                        print(f"Number of iterations: {bisection_iterations}")
                except ValueError as e:
                    print(f"Bisection method failed: {e}")

            else:
                print("Invalid choice. Please try again.")
                break

        print("-" * 50)


if __name__ == "__main__":
    main()
