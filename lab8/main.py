import numpy as np

def Choletsky(A):
    N = A.shape[0]
    L = np.zeros(shape=(N, N))

    L[0][0] = np.sqrt(A[0][0])

    for i in range(1, N):
        L[i][0] = A[i][0] / L[0][0]
    
    for i in range(1, N):
        for j in range(1, i):
            L[i][j] = (1/L[j][j])*(A[i][j] - np.sum(L[i][:j+1]*L[j][:j+1]))
        L[i][i] = np.sqrt(A[i][i] - np.sum(L[i][:i+1]*L[i][:i+1]))
    return L

def LU(A):
    N = A.shape[0]
    L = np.identity(N)
    U = np.zeros(shape=(N, N))

    for i in range(N):
        for j in range(N):
            if i <= j:
                U[i][j] = A[i][j] - np.sum(L[i, :i+1]*U[:i+1, j], axis=-1)
            else:
                L[i][j] = (A[i][j] - np.sum(L[i, :j+1]*U[:j+1, j], axis=-1)) / U[j][j]

    return L, U

def lin_LU(A, b):
    L, U = LU(A)
    Y = []
    print(b)
    y = b[0]
    Y.append(y)
    for i in range(1, len(b)):
        y = b[i] - np.sum(L[i][:i] * Y[:i])
        Y.append(y)

    X = len(b)*[0]
    for i in range(len(b) - 1, -1, -1):
        if len(b) - 1 == i:
            x = Y[i] / U[i][i]
        else: 
            print("U: ", U[i][i+1:])
            print("X: ", X[i+1:])
            x = (Y[i] - np.sum(U[i][i+1:]*X[i+1:])) / U[i][i]
        X[i] = x


    print(X)
    return 

def main():
    N = 5
    a = 1
    b = 10
    A = np.random.uniform(a, b, (N, N))
    A = np.dot(A, A.T)
    print(A)
    L = Choletsky(A)
    newA = L @ L.T
    print()
    print(newA)
    print("err Холецкий: ", np.linalg.norm(A - newA))

    L, U = LU(A)
    newA = L @ U
    print("err LU: ", np.linalg.norm(A - newA))

    x_answ = np.array([1, 1, 1, 1, 1])
    b = A @ x_answ
    lin_LU(A, b)

main()