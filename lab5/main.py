import numpy as np

def f(x, koef):
    sum = 0
    for i in range(len(koef)):
        sum += koef[i] * x**(len(koef)-1-i)
    return sum

delta = 10**(-3)

def BisectionMethod(xn, xk, F, koef): ## xn - начало отрезка xk - конец отрезка
    eps = 10**(-10) # точность поиска
    x_new = 0
    while abs(xn-xk) > eps:
        x_new = (xn + xk) / 2 # середина

        if np.sign(F(xn, koef)) != np.sign(F(x_new, koef)):
            xk = x_new
        else:
            xn = x_new
    return x_new

def searchSegments(l, R, F, koef):
    Segments = []
    nextPoint = l + delta
    while l < R:
        if np.sign(F(l, koef)) != np.sign(F(nextPoint, koef)):
            Segments.append([l, nextPoint])
            l = nextPoint
        if nextPoint > R:
            break
        nextPoint += delta
    return Segments

def getRoots(l, R, F, koef):
    roots = []
    for i in searchSegments(l, R, F, koef):
        answ = BisectionMethod(i[0], i[1], F, koef)
        roots.append(answ)
        #print(i, "  ", answ)
    return roots

def Gershgorin(matrix):
    circles = []
    for i in range(len(matrix)):
        a = matrix[i][i]
        r = sum(map(abs, matrix[i]))
        circles.append([a, r])
    return circles

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

def mulMatrixVecotor(matrix, vector):
    result = [0] * len(matrix)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            result[i] += matrix[i][j] * vector[j]
    return result



def getXVector(ys, q):
    ys = ys[len(ys)-2::-1]
    sumY = []
    for i in range(len(ys[0])):
        yi = []
        s = 0.0
        for j in range(len(ys)):
            s += ys[j][i] * q[j]
        sumY.append(s)
    return sumY



def normolize(vector):
    sum = 0
    for i in vector:
        sum += i**2
    for i in range(len(vector)):
        vector[i] /= (sum**0.5)
    return vector

def method_Krylova(matrix):
    y0 = [0] * len(matrix)
    y0[0] = 1
    ys = []
    ys.append(y0)
    for i in range(len(matrix)):
        y = mulMatrixVecotor(matrix, ys[i])
        ys.append(y)

    system = []
    for i in range(len(matrix)):
        line = []
        for j in range(len(ys) - 2,  -1, -1):
            line.append(ys[j][i])
        system.append(line)

    print(ys)
    equal = ys[-1].copy()

    YS = ys.copy()

    matrix, equal = gauss_method(system, equal)

    n = len(matrix)

    list_X_new = [0] * n
    
    koefs = [1]
    for i in range(n - 1, -1, -1):
        sum_ = equal[i]
        for j in range(i + 1, n):
            sum_ -= matrix[i][j] * list_X_new[j]
        list_X_new[i] = sum_ / matrix[i][i]
       
    print(list_X_new)
    for i in range(n):
        koefs.append(-list_X_new[i])
    print(koefs)
    circles = Gershgorin(system)

    L, R = 0, 0
    for c in circles:
        if c[0] - c[1] < L:
            L = c[0] - c[1]
        if c[0] + c[1] > R:
            R = c[0] + c[1]
    print(L, R)
    roots = getRoots(L, R, f, koefs)
   # print(roots)
    print("----------------------------\n", YS)
    X = []
    for root in roots:
        q = []
        q.append(1)
        for k in range(1, len(koefs) - 1):
            qi = root*q[k-1] + koefs[k]
            q.append(qi)
        
        X.append(normolize(getXVector(YS, q)))
    
    print(X)
    return roots

def main():
    matrix = [[2.2, 1, 0.5, 2], [1, 1.3, 2, 1], [0.5, 2, 0.5, 1.6], [2, 1, 1.6, 2]] 
    method_Krylova(matrix)   
    return

main()