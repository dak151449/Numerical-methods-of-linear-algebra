import numpy as np
import random
import copy
import matplotlib.pyplot as plt

# from ..lab6.lab6 import matrix_random
from Krylov import method_Krylova
from linalg_methods import *

def normMatrix(A):
    s = []
    for i in range(len(A)):
        row_sum = 0
        for j in range(len(A)):
            if i == j:
                continue
            row_sum += abs(A[i][j])
        row_sum = row_sum / A[i][i]
        s.append(row_sum)
    return max(s)

def normaVect(v):
    s = 0
    for i in list(v):
        s += i**2
    return s**(1/2)
    

def norm(x_old, x_new):
    return np.linalg.norm(x_new - x_old)

def gen_symm_matrix(dim , a, b , diag= 10 ) :
    matrix = np.random.uniform( a , b , size =(dim , dim))
    symmetric_matrix = np.triu(matrix) + np.triu(matrix, 1).T

    symmetric_matrix += np.eye(dim)*diag
    return symmetric_matrix

def generSimmetricalMatrix(n, a, b, diag=10):
    matrix = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        # matrix.append([])
        # l = [0 for i in range(n)]
        for j in range(n):
            if i == j:
                matrix[i][j] = random.uniform(a, b) + b*n
            elif matrix[i][j] == 0:
                matrix[i][j] = random.uniform(a, b)
                matrix[j][i] = matrix[i][j]
    return matrix


def OneParameterMethod(A, t, f, eps, eigVal, xAnswer):
    eigvalsP = np.abs(1 - t*np.array(eigVal))
    eigvalsPmax = max(eigvalsP)

    if eigvalsPmax > 1:
        print("eigvalsPmax > 1")
        return np.inf

    A = np.matrix(A)
    G = np.identity(len(A)) - t*A
    f = np.array(f)
    f = f*t
    x_old = np.array([0 for i in range(len(A))])
    # x_new = np.ones((len(A), 1))
    # print(A, '\n',f, '\n', x_old)
    R0 = np.array(xAnswer) - x_old
    R_back = R0
    it = 0
    while True:
        x_new = G@x_old + f
        x_new = np.array(x_new)[0]

        if normaVect(x_old - x_new) < eps:
            break
        x_old = x_new
        it += 1

    

    for i in range(it):
        R_back = R0
        R0 = np.array(G @ R0)[0]
        if np.round(normaVect(R0)**2, 9) > np.round(eigvalsPmax**2 * normaVect(R_back)**2,9) + 1e-4:
            print("Условие сходимости для тау не выполнилось", np.round(np.linalg.norm(R0)**2, 10) -  np.round(eigvalsPmax**2 * np.linalg.norm(R_back)**2, 10), i, it)
            break
        

    print(x_new)
    return it


def lab7():
    random.seed(100)
    eps = 1e-13
    n = 4
    A = generSimmetricalMatrix(n, 1, 2)
    f = [random.uniform(10, 20) for i in range(n)]
    X, _ = Zeidal(A, f, eps)
    print(A, '\n')

    if normMatrix(A) > 1:
        print("normMatrix > 1")
        return

    roots = method_Krylova(A)
    t_opt = 2 / (min(roots) + max(roots))
    
    l = (2 / max(roots)) / 0.001

    t_range = np.arange(0.001, 2 / max(roots), 0.001)
    it_range  = []
    # t_range = [t]
    print("start OneParameterMethod", l)
    index = 0
    for t in t_range:

        it = OneParameterMethod(A, t, f, eps, roots, X)
        it_range.append(it)
        # print(it)
        if index % 10000 == 0:
            print(index)
        index += 1

    

    index = it_range[::-1].index(min(it_range))
    print("fnd: ", t_range[-1 - index])
    print("opt: ", t_opt)
    print("delta: ", abs(t_range[-1 - index] - t_opt))
    print(X)
    print(OneParameterMethod(A, t_opt, f, eps, roots, X), " ", min(it_range))
    plt.plot(t_range, it_range)
    plt.scatter([t_range[-1 - index]], [min(it_range)])
    plt.scatter([t_opt], [OneParameterMethod(A, t_opt, f, eps, roots, X)])
    plt.show()

lab7()


