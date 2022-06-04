import copy
from math import sqrt

IdentityMatrix = lambda size : [[0 if y != x else 1 for y in range(size)] for x in range(size)]

def norm(r):
    sum = 0
    for i in range(len(r)):
        for j in range(len(r[i])):
            #sum += abs(r[i][j])
            sum += r[i][j]**2
    sum = sqrt(sum)
    return sum

def addMatrices(A, B):
    if len(A) != len(B):
        print("not valid matrices to add")
        return None
    C = [[] for x in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[i])):
            C[i].append(A[i][j] + B[i][j])
    return C


def multiplyMatrices(A, B):
    C = [[] for x in range(len(A))]
    for i in range(len(A)):
        if len(A[i]) != len(B):
            print("not valid matrices")
            return None
        for j in range(len(B[0])):
            sum = 0
            for k in range(len(A[i])):
                sum += A[i][k] * B[k][j]
            C[i].append(sum)
    return C


def Gauss_Jordan(X, y):
    A = copy.deepcopy(X)
    b = copy.deepcopy(y)
    for i in range(len(A)):
        pivot = A[i][i]
        n = 0
        while not pivot and n < len(A):  # switch rows if pivot is 0
            A[i], A[i+n] = A[i+n], A[i]
            b[i], b[i + n] = b[i + n], b[i]
            pivot = A[i][i]
            n += 1

        #b[i] /= pivot
        if pivot != 1:
            for j in range(i, len(A[i])): # divide row by pivot
                A[i][j] /= pivot

            for j in range(0, len(b[i])): # divide row by pivot
                b[i][j] /= pivot

        for k in range(i, len(A)):
            if k != i:
                factor = A[k][i]
                #b[k][i] = b[k][i] - factor * b[i][i]
                if factor:
                    for j in range(i,len(A[k])):
                        A[k][j] = A[k][j] - factor * A[i][j]
                    for j in range(0,len(b[k])):
                        b[k][j] = b[k][j] - factor * b[i][j]
    return A, b


def inverseMatrix(A):
    I = IdentityMatrix(len(A))
    I, A1 = Gauss_Jordan(A,I)
    return A1


def multiplyByNum(A, n):
    X = [[] for x in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[i])):
            if A[i][j]:
                X[i].append(A[i][j] * n)
            else:
                X[i].append(0)
    return X


def Jacobi(L, U, D, b):
    r = [[1] for x in range(N)] # r0
    LU = addMatrices(L,U)
    M = addMatrices(LU, D)
    minusB = multiplyByNum(b,-1) # -b
    D1 = inverseMatrix(D) # D(-1)
    D1_b = multiplyMatrices(D1,b) # D(-1)b
    D1_LU = multiplyMatrices( multiplyByNum(D1, -1), LU) # -D(-1)(L+U)

    start = time()
    r = addMatrices(multiplyMatrices(D1_LU ,r), D1_b) # r = -D(-1)(L+U)r + D(-1)b
    res = addMatrices(multiplyMatrices(M,r), minusB) # res = Mr - b
    iteration = 1
    while norm(res) > RES_TRESHHOLD:
        r = addMatrices(multiplyMatrices(D1_LU ,r), D1_b)
        res = addMatrices(multiplyMatrices(M, r), minusB)
        iteration += 1
    print(f'Jacobi iterations: {iteration}')
    end = time()
    print(f'time = {end-start}')
    return r


def Gauss_Siedel(L, U, D, b):
    r = [[1] for x in range(N)] # r0
    d = D
    M = addMatrices(addMatrices(L,U), d)
    minusB = multiplyByNum(b,-1) # -b
    DL = addMatrices(D,L)
    DL1 = inverseMatrix(DL)
    DL1b = multiplyMatrices(DL1,b)
    DL1U = multiplyMatrices(multiplyByNum(DL1,-1), U)

    start = time()
    r = addMatrices(multiplyMatrices(DL1U, r), DL1b)
    res = addMatrices(multiplyMatrices(M,r), minusB)
    iteration = 1
    while norm(res) > RES_TRESHHOLD:
        r = addMatrices(multiplyMatrices(DL1U, r), DL1b)
        res = addMatrices(multiplyMatrices(M, r), minusB)
        #print(f'{iteration}: {norm(res)}')
        iteration += 1
    print(f'Gauss-Siedel iterations: {iteration}')
    end = time()
    print(f'time = {end-start}')
    return r
