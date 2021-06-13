def precision(num):
    epsilon = 1.0
    while (1.0 + 0.5 * epsilon) != 1:
        epsilon = 0.5 * epsilon
    if abs(num) <= epsilon:
        return 0
    else:
        return num


def vector_mul(vector, matrix):
    result = [[0] * len(matrix)]
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            result[0][i] += matrix[i][j] * vector[j]
    return result


def MatrixAddition(matrix_1, matrix_2):
    return [[matrix_1[i][k] + matrix_2[i][k] for k in range(len(matrix_1))] for i in range(len(matrix_1))]


def matrixMultiply(matrix_1, matrix_2):
    result = [[0] * len(matrix_2) for i in range(len(matrix_1))]  # Initializing the result matrix (with zeros)
    for i in range(len(matrix_1)):
        for j in range(len(matrix_2)):
            for k in range(len(matrix_2)):
                result[i][j] += precision(matrix_1[i][k] * matrix_2[k][j])
    return result


def createMatrix(matrix, size, val):
    if val == 0:
        matrix = [[0 for i in range(size)] for j in range(size)]
    elif val == 1:
        matrix = [[int(i == j) for i in range(size)] for j in range(size)]
    return matrix


def calcDet(matrix):
    size = len(matrix)
    det = 0
    if size == 2:
        det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        return det
    minor = createMatrix(matrix, size - 1, 0)
    for k in range(len(matrix)):
        i, j = 0, 0
        while i < size:
            if i != k:
                minor[j] = matrix[i][1:]
                j += 1
            i += 1
        det += matrix[k][0] * ((-1) ** k) * calcDet(minor)
    return det


def invertMatrix(matrix):
    determinant = calcDet(matrix)
    inverse = []
    for i in range(len(matrix)):
        inverseRow = []
        for j in range(len(matrix)):
            minor = [row[:j] + row[j + 1:] for row in (matrix[:i] + matrix[i + 1:])]
            inverseRow.append(((-1) ** (i + j)) * calcDet(minor))
        inverse.append(inverseRow)
    inverse = list(map(list, zip(*inverse)))
    for i in range(len(inverse)):
        for j in range(len(inverse)):
            inverse[i][j] = inverse[i][j] / determinant
    return inverse


def D(matrix):
    D = []
    D = createMatrix(D, len(matrix), 0)
    for i in range(len(matrix)):
        D[i][i] = matrix[i][i]
    return D


def L(matrix):
    L = []
    L = createMatrix(L, len(matrix), 0)
    for i in range(len(matrix)):
        for n in range(i):
            L[i][n] = matrix[i][n]
    return L


def U(matrix):
    U = []
    U = createMatrix(U, len(matrix), 0)
    for i in range(len(matrix)):
        for m in range(i + 1, len(matrix)):
            U[i][m] = matrix[i][m]
    return U


def gauss_G(D, L, U):
    G = matrixMultiply(invertMatrix(MatrixAddition(L, D)), U)
    for i in range(len(G)):
        for j in range(len(G)):
            G[i][j] = -G[i][j]
    return G


def norma(G):
    m_sum = 0
    for i in range(len(G)):
        i_sum = 0
        for j in range(len(G)):
            i_sum += abs(G[i][j])
            if i_sum > m_sum:
                m_sum = i_sum
    return m_sum


def hasDominantDiagonal(matrix):
    for i in range(len(matrix)):
        pivot, sum = abs(matrix[i][i]), 0
        for k in range(len(matrix)):
            sum = sum + abs(matrix[i][k])
        if sum - pivot > pivot:
            return False
    return True


def gaussSeidelMethod(matrix, b):
    Xr, Yr, Zr, condition, count, epsilon = 0, 0, 0, 1, 0, 0.00001

    while condition > epsilon and count < 100:
        count += 1
        Xr_1 = Xr
        Xr = (b[0] - matrix[0][1] * Yr - matrix[0][2] * Zr) / matrix[0][0]
        Yr = (b[1] - matrix[1][0] * Xr - matrix[1][2] * Zr) / matrix[1][1]
        Zr = (b[2] - matrix[2][1] * Yr - matrix[2][0] * Xr) / matrix[2][2]
        condition = abs(Xr - Xr_1)
        Xr, Yr, Zr = round(Xr, 6), round(Yr, 6), round(Zr, 6)
        print("{0} X = {1}, Y = {2}, Z = {3}".format(count, Xr, Yr, Zr))
    return count, Xr, Yr, Zr


def driver(matrix, b):
    if hasDominantDiagonal(matrix):
        print("\nDominant Diagonal Matrix: Yes")
    else:
        print("\nDominant Diagonal Matrix: No")

    print('\n*****Gauss Elimination Method*****')
    print('X = {0}, Y = {1}, Z = {2}'.format(vector_mul(b, invertMatrix(matrix))[0][0],vector_mul(b, invertMatrix(matrix))[0][1],vector_mul(b, invertMatrix(matrix))[0][2]))

    print("\n****Gauss Seidel Method****")
    gaussSeidelMethod(matrix, b)
    if norma(gauss_G(D(matrix), L(matrix), U(matrix))) < 1:
        print("Converge !")
    else:
        print("Not converge...")


A = [[2, 1, 0], [3, -1, 0], [1, 4, -2]]
B = [-3, 1, -5]

driver(A, B)
