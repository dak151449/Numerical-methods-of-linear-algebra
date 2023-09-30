import numpy as np
import random
import copy
import matplotlib.pyplot as plt

def norm(a):
  return np.sqrt(sum(x ** 2 for x in a))

def mul_matrix_vec(a, b):
  N = len(a)
  M = len(b)

  c = [0] * N
  for i in range(N):
    for j in range(M):
      c[i] += a[i][j] * b[j]
  return c

def copyM(M):
  n = len(M)
  M1 = [[0] * n for i in range(n)]
  for i in range(n):
    for j in range(n):
      M1[i][j] = M[i][j]
  return M1


def GaussRows(a2, f2):
  a = copyM(a2)
  new_f = copy.copy(f2)
  n = len(a)
  x_order = [i for i in range(n)]

  for k in range(0, n - 1):

    # findMax:
    #print(a)
    max = abs(a[k][k])
    maxind = k
    for i in range(k, n):
      if abs(a[k][i]) > max:
        max = abs(a[k][i])
        maxind = i

    for i in range(n):
      a[i][k], a[i][maxind] = a[i][maxind], a[i][k]

    x_order[k], x_order[maxind] = x_order[maxind], x_order[k]
    #print(x_order)
    #print(a)


    for i in range(k+1, n):
      l = -a[i][k] / a[k][k]
      new_f[i] = new_f[i] + l * new_f[k]
      a[i][k] = 0
      for j in range (k+1, n):
        a[i][j] = a[i][j] + l*a[k][j]

  new_x = [0] * n
  for i in range(n-1, -1, -1):
    sum = new_f[i]
    for j in range(i+1, n):
      sum -= a[i][j] * new_x[j]
    new_x[i] = sum / a[i][i]

  x = [0] * n
  for i in range(len(x_order)):
    x[x_order[i]] = new_x[i]
  return x


def GaussColon(a2, f2):
  a = copyM(a2)
  new_f = copy.copy(f2)
  n = len(a)
  x_order = [i for i in range(n)]

  for k in range(0, n - 1):

    max = abs(a[k][k])
    maxind = k

    for i in range(k, n):
      if abs(a[i][k]) > max:
        max = abs(a[i][k])
        maxind = i

    for i in range(n):
      a[k][i], a[maxind][i] = a[maxind][i], a[k][i]

    new_f[k], new_f[maxind] = new_f[maxind], new_f[k]



    for i in range(k+1, n):
      l = -a[i][k] / a[k][k]
      new_f[i] = new_f[i] + l * new_f[k]
      a[i][k] = 0
      for j in range (k+1, n):
        a[i][j] = a[i][j] + l*a[k][j]

  new_x = [0] * n
  for i in range(n-1, -1, -1):
    sum = new_f[i]
    for j in range(i+1, n):
      sum -= a[i][j] * new_x[j]
    new_x[i] = sum / a[i][i]

  return new_x


def GaussComb(a2, f2):
  a = copyM(a2)
  new_f = copy.copy(f2)
  n = len(a)
  x_order = [i for i in range(n)]

  for k in range(0, n - 1):

    # findMax:
    #print(a)
    max = abs(a[k][k])
    maxind = k
    swap = 'colon'

    for i in range(k, n):
      if abs(a[i][k]) > max:
        max = abs(a[i][k])
        maxind = i
    for i in range(k, n):
      if abs(a[k][i]) > max:
        max = abs(a[k][i])
        maxind = i
        swap = 'row'

    if swap == 'row':
      for i in range(n):
        a[i][k], a[i][maxind] = a[i][maxind], a[i][k]
      x_order[k], x_order[maxind] = x_order[maxind], x_order[k]
    else:
      for i in range(n):
        a[k][i], a[maxind][i] = a[maxind][i], a[k][i]
      new_f[k], new_f[maxind] = new_f[maxind], new_f[k]

    '''print(x_order)
    print(a)'''


    for i in range(k+1, n):
      l = -a[i][k] / a[k][k]
      new_f[i] = new_f[i] + l * new_f[k]
      a[i][k] = 0
      for j in range (k+1, n):
        a[i][j] = a[i][j] + l*a[k][j]

  new_x = [0] * n
  for i in range(n-1, -1, -1):
    sum = new_f[i]
    for j in range(i+1, n):
      sum -= a[i][j] * new_x[j]
    new_x[i] = sum / a[i][i]

  x = [0] * n
  for i in range(len(x_order)):
    x[x_order[i]] = new_x[i]
  return x


def count_err(new_x, x):
  n = len(new_x)
  dx = [0]*n
  for i in range(n):
    dx[i] = abs(new_x[i] - x[i])
  return norm(dx)

def matrix_random(N, A, B, diag):
    Matrix = [[random.randint(A, B) for j in range(N)] for i in range(N)]

    for i in range(N):
        Matrix[i][i] = sum(list(map(abs, Matrix[i]))) + diag - abs(Matrix[i][i])

    return Matrix


def vector_random(N, A, B):
  v = [0] * N
  for i in range(N):
    r = random.randint(A, B)
    v[i] = r
  return v


def main():
  n = 10
  print("\nrandom: ", n)

  x_graph = [ i for i in range(-40, 21)]# if i != 0]
  y1 = []
  y2 = []
  y3 = []
  y4 = []

  for i in range(-40, 21):
    ##if i == 0: continue
    M = matrix_random(n, -10, 10, i)

    x = vector_random(n, -10, 10)
    f = mul_matrix_vec(M, x)


    x_rows = GaussRows(M, f)
    x_colon = GaussColon(M, f)
    x_comb = GaussComb(M, f)
    x_lib = np.linalg.solve(M, f)
    y1.append(count_err(x_rows, x))
    y2.append(count_err(x_colon, x))
    y3.append(count_err(x_comb, x))
    y4.append(count_err(x_lib, x))

  plt.plot(x_graph, y1, label='rows')
  plt.plot(x_graph, y2, label='colons')
  plt.plot(x_graph, y3, label='comb')
  plt.plot(x_graph, y4, label='lib')
  plt.grid()
  plt.legend()
  plt.show()

main()
