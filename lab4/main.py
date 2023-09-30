import random
import numpy as np
import matplotlib.pyplot as plt

def funcShow(f):
    x = np.linspace(-10, 10, 2000)
    plt.plot(x, f(x))
    plt.plot(x, 0*x)
    plt.show()

koef = [1, 0, -8, 1]

def f(x):
    sum = 0
    for i in range(len(koef)):
        #print(koef[i])
        sum += koef[i] * x**(len(koef)-1-i)
    #print()
    return sum

# def f(x):
#     return x**3 - 8*x + 1

delta = 10**(-3)



def BisectionMethod(xn, xk, F): ## xn - начало отрезка xk - конец отрезка
    eps = 10**(-10) # точность поиска
    x_new = 0
    while abs(xn-xk) > eps:
        x_new = (xn + xk) / 2 # середина

        if np.sign(F(xn)) != np.sign(F(x_new)):
            xk = x_new
        else:
            xn = x_new
    return x_new

def searchSegments(l, R, F):
    Segments = []
    nextPoint = l + delta
    while l < R:
        if np.sign(F(l)) != np.sign(F(nextPoint)):
            Segments.append([l, nextPoint])
            l = nextPoint
        if nextPoint > R:
            break
        nextPoint += delta
    return Segments

def getRoots(l, R, F):
    roots = []
    for i in searchSegments(l, R, F):
        answ = BisectionMethod(i[0], i[1], F)
        roots.append(answ)
        print(i, "  ", answ)
    return roots

def eigenValue(matrix):
    values = []
    for i in range(len(matrix)):
        values.append(matrix[i][i])
    return values


def main():
    getRoots(-100, 100, f)

    funcShow(f)
    
    return

main()
