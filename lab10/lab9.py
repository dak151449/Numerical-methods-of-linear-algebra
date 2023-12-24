import numpy as np

lmbda = []
p = []

def f(x):
  n = len(p) - 1
  sum = 0
  for i in range (len(p)):
    sum += - (-1)**n * p[i] * x**(n-i)
  return sum


def find0(a, b):
  eps = 10e-5
  if f(a) < eps:
    lmbda.append(a)
    return
  if f(b) < eps:
    lmbda.append(b)
    return
  x1 = (a + b) / 2
  if f(x1)*f(a) <= 0:
    find0(a, x1)
  if f(x1)*f(b) <= 0:
    find0(x1, b)

def resolveEq(A, B):
  n = 10000
  delta = (B-A) / n
  a = A
  for i in range(n):
    b = a + delta
    if f(a)*f(b) <= 0:
      find0(a, b)
    a = b

def getGershgorin(M):
  n = len(M)
  left = 1000
  right = -1000
  for i in range(n):
    r = 0;
    z = M[i][i]
    for j in range(n):
      if j != i:
        r += abs(M[i][j])
    left = (z - r) if left > z-r else left
    right = (z + r) if right < z+r else right

  return left, right

def Danilevski(M):
  n = len(M)
  D = [row[:] for row in M]
  B = []
  for k in range (n-2, -1, -1):
    B1 = [[0 if i!=j else 1 for j in range(n)] for i in range(n)]
    for i in range(n):
      if i!= k:
        B1[k][i] = - D[k+1][i] / D[k+1][k]
      else:
        B1[k][i] = 1 / D[k+1][k]

    B1_inv = np.linalg.inv(B1)
    D = np.dot(np.dot(B1_inv, D), B1)
    if k == n-2:
      B = [row[:] for row in B1]
    else:
      B = np.dot(B, B1)

  for d in D[0]:
    p.append(d)
  return B

def findVectors(B):
  n = len(B)
  xs = []
  for l in lmbda:
    y = [l**i for i in range(n-1, -1, -1)]
    x = np.dot(B, y)
    xs.append(x)
  return xs

def ortonorm(xs):
  new_x = []
  for x in xs:
    x_norm = np.linalg.norm(x)
    new_x.append([x[i]/x_norm for i in range(len(x))])
  return new_x

def getDanilevski(A):
  lmbda.clear()
  p.clear()
  p.append(-1)
  B = Danilevski(A)
  l, r = getGershgorin(A)
  resolveEq(l, r)
  x = findVectors(B)
  x = ortonorm(x)
  return np.array(lmbda), np.array(x)



def svd(A):
    AtA = np.dot(A.T, A)
    eigenvalues, V = np.linalg.eig(AtA)

    S = np.diag(np.sqrt(eigenvalues))
    Sinv = [1/s for s in np.sqrt(eigenvalues)]
    Sinv = np.diag(Sinv)

    U = np.dot(np.dot(A, V), Sinv)
    return U, S, V

# A = np.array([[3, 1, 0], [1, 2, 2], [0, 1, 1]])
N = 30
a = 1
b = 10
A = np.random.uniform(a, b, (N, N))

# print("Матрица А:")
# print(A)
# print()

U, S, V = svd(A)
# print("Левые сингулярные векторы (U):")
# print(U)
# print("Сингулярные значения (S):")
# print(S)
# print("Правые сингулярные векторы (V):")
# print(V)
# print("----------------\n", len(U), len(U[0]), len(V), len(V[0]), len(np.eye(N))
print(np.linalg.norm((U@V.T)@((U@V.T).T) - np.eye(N)))
# print(U@V.T)
res = np.dot(np.dot(U, S), V.T)
print("Проверка перемножением:")
# print(res)
print("mysvd norm dif: ", np.linalg.norm(A - res))

print()
print("numpy:")

# numpy
U, S, Vt = np.linalg.svd(A)
S = np.diag(S)

# print("Левые сингулярные векторы (U):\n", U)
# print("Сингулярные значения (S):\n", S)
# print("Правые сингулярные векторы (VT):\n", Vt)

print(np.linalg.norm((U@Vt)@((U@Vt).T) - np.eye(N)))
numpyA = np.dot(np.dot(U, S), Vt)
print("numpy norm dif: ", np.linalg.norm(A - numpyA))