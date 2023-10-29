import math
import random

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
            equal[k] += c * equal[i]

    return matrix, equal

def norm(X_old, X_new):
    return math.sqrt(sum([(X_new[i]-X_old[i])**2 for i in range(len(X_old))]))

def matrix_random(N, A, B, diag):
    answer = [random.randint(A, B) for j in range(N)]
    Matrix = [[random.randint(A, B) for j in range(N)] for i in range(N)]

    for i in range(N):
        Matrix[i][i] = sum(list(map(abs, Matrix[i]))) + diag - abs(Matrix[i][i])

    return Matrix, answer

def GetDiagonal(A, answer, X_old):
    X_new = []
    for i in range(len(A)):
        sum = answer[i]
        for j in range(len(A[0])):
            if (i == j):
                continue   
            sum -= A[i][j]*X_old[j]

        sum = sum / A[i][i]
        X_new.append(sum)
    return X_new

def GetDiagonalZeidal(A, answer, X_old):
    X_new = []
    for i in range(len(A)):
        sum = answer[i]
        for j in range(len(A[0])):
            if (i == j):
                continue 
            if (j < len(X_new)):
                sum -= A[i][j]*X_new[j]
            else:
                sum -= A[i][j]*X_old[j]

        sum = sum / A[i][i]
        X_new.append(sum)
    return X_new

def Yakoby(matrix, answer, eps):
    X = [0 for i in range(len(matrix))]
    iter_count = 0
    while True:
        X_new = GetDiagonal(matrix, answer[::], X)
        iter_count += 1
        if norm(X, X_new) < eps:
            X = X_new
            break
        X = X_new

    print("Yakoby: ", X)
    print("Iter_count: ", iter_count, "Eps: ", eps)
    return X, iter_count

def Zeidal(matrix, answer, eps):
    X = [0 for i in range(len(matrix))]
    iter_count = 0
    while True:
        X_new = GetDiagonalZeidal(matrix, answer[::], X)
        iter_count += 1
        if norm(X, X_new) < eps:
            X = X_new
            break
        X = X_new

    print("Zeidal: ", X)
    print("Iter_count: ", iter_count, "Eps: ", eps)
    return X, iter_count

def lab6():
    eps = 1e-16
    matrix, answer  = matrix_random(8, -10, 10, 1) ##[[1,200,3], [4,500,6], [7,800,9]]
    
    x_Yakoby, i_Y = Yakoby(matrix, answer, eps)
    x_Zeidal, i_Z = Zeidal(matrix, answer, eps)
    
    matrix, equal = gauss_method(matrix, answer)

    n = len(matrix)

    list_X_new = [0] * n
    
    koefs = [1]
    for i in range(n - 1, -1, -1):
        sum_ = equal[i]
        for j in range(i + 1, n):
            sum_ -= matrix[i][j] * list_X_new[j]
        list_X_new[i] = sum_ / matrix[i][i]
       
    print("Gauss: ", list_X_new, "\nYakoby_norm", norm(x_Yakoby, list_X_new), "\nZeidal_norm", norm(x_Zeidal, list_X_new))
    print("Iter_Y / Iter_Z: ", i_Y / i_Z)
    return

lab6()