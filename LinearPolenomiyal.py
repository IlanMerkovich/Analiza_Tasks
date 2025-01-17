from MatrixUtility import *

def linearInterpolation(table_points, point):
    """
    Linear interpolation and extrapolation function.
    :param table_points: List of data points
    :param point: Value to interpolate or extrapolate
    :return: Interpolation or extrapolation result
    """
    p = []
    result = 0
    flag = 1
    for i in range(len(table_points)):
        p.append(table_points[i][0])
    for i in range(len(p) - 1):
        if i <= point <= i + 1:
            x1 = table_points[i][0]
            x2 = table_points[i + 1][0]
            y1 = table_points[i][1]
            y2 = table_points[i + 1][1]
            result = (((y1 - y2) / (x1 - x2)) * point) + ((y2 * x1) - (y1 * x2)) / (x1 - x2)
            print("\nThe approximation (interpolation) of the point ", point, " is: ", round(result, 4))
            flag = 0
    if flag:
        x1 = table_points[0][0]
        x2 = table_points[1][0]
        y1 = table_points[0][1]
        y2 = table_points[1][1]
        m = (y1 - y2) / (x1 - x2)
        result = y1 + m * (point - x1)
        print("\nThe approximation (extrapolation) of the point ", point, " is: ", round(result, 4))
    return result


def GaussJordanElimination(matrix, vector):
    """
    Function for solving a linear equation using Gauss's elimination method
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=A(-1)b
    """
    matrix, vector = RowXchange(matrix, vector)
    invert = InverseMatrix(matrix, vector)
    return MulMatrixVector(invert, vector)


def UMatrix(matrix, vector):
    """
    :param matrix: Matrix nxn
    :return: Disassembly into a U matrix
    """
    U = MakeIMatrix(len(matrix), len(matrix))
    for i in range(len(matrix[0])):
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)
    U = MultiplyMatrix(U, matrix)
    return U


def LMatrix(matrix, vector):
    """
    :param matrix: Matrix nxn
    :return: Disassembly into an L matrix
    """
    L = MakeIMatrix(len(matrix), len(matrix))
    for i in range(len(matrix[0])):
        matrix, vector = RowXchageZero(matrix, vector)
        for j in range(i + 1, len(matrix)):
            elementary = MakeIMatrix(len(matrix[0]), len(matrix))
            elementary[j][i] = -(matrix[j][i]) / matrix[i][i]
            L[j][i] = (matrix[j][i]) / matrix[i][i]
            matrix = MultiplyMatrix(elementary, matrix)
    return L


def SolveLU(matrix, vector):
    """
    Function for solving a linear equation by LU decomposition
    :param matrix: Matrix nxn
    :param vector: Vector n
    :return: Solve Ax=b -> x=U(-1)L(-1)b
    """
    matrixU = UMatrix(matrix, vector)
    matrixL = LMatrix(matrix, vector)
    return MultiplyMatrix(InverseMatrix(matrixU), MultiplyMatrix(InverseMatrix(matrixL), vector))


def solveMatrix(matrixA, vectorb):
    detA = Determinant(matrixA, 1)
    print("\nDET(A) =", detA)

    if detA != 0:
        print("\nNon-Singular Matrix - Perform GaussJordanElimination")
        result = GaussJordanElimination(matrixA, vectorb)
        print(result)
        return result
    else:
        print("Singular Matrix - Perform LU Decomposition\n")
        print("Matrix U:\n", UMatrix(matrixA, vectorb))
        print("\nMatrix L:\n", LMatrix(matrixA, vectorb))
        print("\nMatrix A=LU:\n", MultiplyMatrix(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb)))
        return MultiplyMatrix(LMatrix(matrixA, vectorb), UMatrix(matrixA, vectorb))


def polynomialInterpolation(table_points, x):
    """
    Polynomial interpolation function
    :param table_points: List of data points
    :param x: Value to interpolate
    :return: Interpolation result
    """
    matrix = [[point[0] ** i for i in range(len(table_points))] for point in table_points]
    b = [[point[1]] for point in table_points]

    print("The matrix obtained from the points:\n", matrix)
    print("\nb vector:\n", b)

    matrixSol = solveMatrix(matrix, b)
    result = sum([matrixSol[i][0] * (x ** i) for i in range(len(matrixSol))])

    print("\nThe polynomial:")
    print('P(X) = ' + ' + '.join([f"({matrixSol[i][0]}) * x^{i}" for i in range(len(matrixSol))]))
    print(f"\nThe result of P(X={x}) is:", result)
    return result


if __name__ == '__main__':
    table_points = [(0, 0), (1, 0.8415), (2, 0.9093), (3, 0.1411), (4, -0.7568), (5, -0.9589), (6, -0.2794)]
    x = 1.28

    print("----------------- Interpolation & Extrapolation Methods -----------------\n")
    print("Table Points:")
    for point in table_points:
        print(f"({point[0]}, {point[1]})")
    print("\nFinding an approximation to the point:", x)

    # Linear Interpolation
    print("\n--- Linear Interpolation ---")
    linear_result = linearInterpolation(table_points, x)
    print(f"Linear Interpolation Result: {linear_result}")

    # Polynomial Interpolation
    print("\n--- Polynomial Interpolation ---")
    polynomial_result = polynomialInterpolation(table_points, x)
    print(f"Polynomial Interpolation Result: {polynomial_result}")

    print("\n---------------------------------------------------------------------------\n")
