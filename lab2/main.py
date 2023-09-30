import random
import numpy as np

def mul_matrix(A, B):
    n = len(A)
    m = len(B[0])
    C = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for k in range(len(B)):
                C[i][j] += A[i][k] * B[k][j]
    return C


def gener_matrix(maxx, minn, n, m):
    matrix = n * [0]
    for i in range(n):
        l = []
        for j in range(m):
            l.append(random.uniform(float(minn), float(maxx)))
        matrix[i] = l
    return matrix

def add_diag(matrix, n):
    for i in range(len(matrix)):
        matrix[i][i] += n
    return matrix


def gauss_method(matrix, equal):
    n = len(matrix)
    for i in range(n):
        maxElem = abs(matrix[i][i])
        maxRow = i
        for k in range(i + 1, n):
            if abs(matrix[k][i]) > maxElem:
                maxElem = abs(matrix[k][i])
                maxRow = k

        matrix[i], matrix[maxRow] = matrix[maxRow], matrix[i]
        equal[i], equal[maxRow] = equal[maxRow], equal[i]

        for k in range(i + 1, n):
            c = -matrix[k][i] / matrix[i][i]
            for j in range(i, n):
                if i == j:
                    matrix[k][j] = 0
                else:
                    matrix[k][j] += c * matrix[i][j]
            equal[k][0] += c * equal[i][0]

    return matrix, equal

def gener_3_matrix(n):
    M = []
    for i in range(n):
        M.append([])
        for j in range(n):
            if abs(i - j) <= 1:
                M[i].append(random.uniform(-100, 100))
            else:
                M[i].append(0)
    return M


def progonka(matrix, D, index_line, out):
    if index_line == 0:
        out = [[matrix[0][0], matrix[0][1]/matrix[0][0], D[0][0]/matrix[0][0]]]
        out = progonka(matrix, D, index_line + 1, out)
    elif index_line == len(matrix) - 1:
        yn = matrix[index_line][len(matrix[index_line]) - 1] + matrix[index_line][-2]*out[-1][1]
        bn = (D[index_line][0] -  matrix[index_line][-2]*out[-1][2]) / yn
        out.append([yn, 0, bn])
        return out
    else:
        yi = matrix[index_line][index_line] + matrix[index_line][index_line - 1]*out[-1][1]
        ai = -matrix[index_line][index_line + 1] / yi
        bi = (D[index_line][0] - matrix[index_line][index_line - 1]*out[-1][2]) / yi
        out.append([yi, ai, bi])
        out = progonka(matrix, D, index_line + 1, out)
    return out

def get_ansver_prog(matrix, D):
    out = progonka(matrix, D, 0, [])
    X = [0 for _ in range(len(matrix))]
    out[-1]
    X[-1] = out[-1][2]
    for i in range(len(matrix) - 2, -1, -1):
        X[i] = out[i][1] * X[i + 1] + out[i][2]
    return X

def tridiagonal_solver(A, b):
        A = np.matrix(A)
        n = len(A)
        P = np.zeros(n)
        Q = np.zeros(n)

        P[0] = -A[0, 1] / A[0, 0]
        Q[0] = b[0][0] / A[0, 0]

        for i in range(1, n-1):
            denominator = A[i, i] + A[i, i - 1] * P[i - 1]
            P[i] = -A[i, i + 1] / denominator
            Q[i] = (b[i][0] - A[i, i - 1] * Q[i - 1]) / denominator
        
        x = np.zeros(n)
        x[n - 1] = Q[n - 1]

        for i in range(n - 2, -1, -1):
            x[i] = P[i] * x[i + 1] + Q[i]

        return x


def evklid_norm(list_X, list_Y):
    norm = 0
    for i in range(len(list_X)):
        norm += (list_X[i] - list_Y[i]) ** 2
    return norm ** 0.5

def main():
    n = int(input())
    m = n

    

    matrix = gener_matrix(-100, 100, n, m)
    matrix = gener_3_matrix(n)
    print(matrix)
    list_X = []
    for i in range(n):
        list_X.append([random.uniform(-100,  100)])
    #list_X = [[1/3], [1/3]]
    #print(matrix)
    #matrix = [list(map(int, input().split())) for _ in range(n)]
    #print(matrix)
    equal1 = mul_matrix(matrix, list_X)

    #print(matrix, equal)
    matrix, equal = gauss_method(matrix, equal1)

    list_X_new = [0] * n
    

    for i in range(n - 1, -1, -1):
        sum_ = equal[i][0]
        for j in range(i + 1, n):
            sum_ -= matrix[i][j] * list_X_new[j]
        list_X_new[i] = sum_ / matrix[i][i]

    #print(list_X_new)
    l = []
    for i in range(n):
        l.append(list_X[i][0])
    print("Гаусс: ", evklid_norm(l, list_X_new))

    #answer = get_ansver_prog(matrix, equal1)
    answer = tridiagonal_solver(matrix, equal1)
    print(l, '\n', answer, '\n', list_X_new)
    print("Прогонка: ", evklid_norm(l, answer)) #evklid_norm(l, answer))

    l = []
    for i in range(n):
        l.append(equal[i][0])
    out = np.linalg.solve(matrix, l)
    #print(out)

    print("np.linalg.solve: ", *evklid_norm(out, list_X))

    
main()