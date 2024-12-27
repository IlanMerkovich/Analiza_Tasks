def is_diagonally_dominant(matrix):
    n = len(matrix)
    for i in range(n):
        row_sum = sum(abs(matrix[i][j]) for j in range(n) if i != j)
        if abs(matrix[i][i]) < row_sum:
            return False
    return True


def make_diagonally_dominant(matrix, vector):
    n = len(matrix)
    indices = list(range(n))
    for i in range(n):
        max_row = max(indices, key=lambda k: abs(matrix[k][i]))
        if abs(matrix[max_row][i]) < sum(abs(matrix[max_row][j]) for j in range(n) if j != i):
            continue
        matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
        vector[i], vector[max_row] = vector[max_row], vector[i]
    return matrix, vector


def process_matrix(matrix, vector):
    if is_diagonally_dominant(matrix):
        print("The matrix is already diagonally dominant.")
        return matrix, vector

    print("The matrix is not diagonally dominant. Attempting to rearrange...")
    matrix, vector = make_diagonally_dominant(matrix, vector)

    if is_diagonally_dominant(matrix):
        print("The matrix has been rearranged and is now diagonally dominant.")
        return matrix, vector
    else:
        print("The matrix cannot be made diagonally dominant.")
        return None, None


def jacobi_method(A, B, tolerance=0.00001, max_iterations=100):
    n = len(A)
    x = [0] * n
    x_new = x[:]
    results = []

    for iteration in range(max_iterations):
        for i in range(n):
            s = sum(A[i][j] * x[j] for j in range(n) if j != i)
            x_new[i] = (B[i][0] - s) / A[i][i]
        results.append(x_new[:])
        if all(abs(x_new[i] - x[i]) < tolerance for i in range(n)):
            return results, iteration + 1
        x = x_new[:]

    return results, max_iterations


def gauss_seidel_method(A, B, tolerance=0.00001, max_iterations=100):
    n = len(A)
    x = [0] * n
    results = []

    for iteration in range(max_iterations):
        x_new = x[:]
        for i in range(n):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (B[i][0] - s1 - s2) / A[i][i]
        results.append(x_new[:])
        if all(abs(x_new[i] - x[i]) < tolerance for i in range(n)):
            return results, iteration + 1
        x = x_new[:]

    return results, max_iterations


def display_results(method_name, results, iterations):
    print(f"\n{method_name} Results:")
    for i, result in enumerate(results, 1):
        print(f"Iteration {i}: {result}")
    print(f"Total iterations: {iterations}")


def main():
    matrixA = [[4, 2, 0], [2, 10, 4], [0, 4, 5]]
    vectorB = [[2], [6], [5]]
    matrixA_processed, vectorB_processed = process_matrix(matrixA, vectorB)
    if not matrixA_processed:
        print("Exiting: The matrix cannot be made diagonally dominant.")
        return
    while True:
        print("\nChoose a method:")
        print("1. Jacobi Method")
        print("2. Gauss-Seidel Method")
        print("3. Exit")
        choice = input("Enter your choice (1/2/3): ").strip()

        if choice == "1":
            results, iterations = jacobi_method(matrixA_processed, vectorB_processed)
            display_results("Jacobi Method", results, iterations)
        elif choice == "2":
            results, iterations = gauss_seidel_method(matrixA_processed, vectorB_processed)
            display_results("Gauss-Seidel Method", results, iterations)

if __name__ == "__main__":
    main()
