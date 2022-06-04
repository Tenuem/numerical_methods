from matrix_operations import multiplyMatrices, IdentityMatrix, addMatrices, multiplyByNum, inverseMatrix, norm, Gauss_Jordan
from math import sin, pow, sqrt
from time import time
#indeks 184640
cd = 40
e = 6
f = 4

a1, a2, a3 = 5+e, -1, -1
N = 900 + cd
RES_TRESHHOLD = 10 ** (-9)

import copy
def Jacobi_method(A, b):
    r = [[0] for x in range(len(A))]  # r0
    X = [[1] for x in range(len(A))]  # r(i+1)
    res = addMatrices(multiplyMatrices(A, r), multiplyByNum(b, -1))
    iteration = 0
    norm_res = norm(res)
    while norm_res > RES_TRESHHOLD:
        if norm_res == float('inf'):
            return "no solution found"
        for i in range(len(A)):
            x = 0
            for j in range(len(A[i])):
                if j != i:
                    x += A[i][j]*r[j][0]
            X[i][0] = (1/A[i][i])*(b[i][0] - x)
        r = copy.deepcopy(X)
        res = addMatrices(multiplyMatrices(A, r), multiplyByNum(b, -1))
        norm_res = norm(res)
        iteration += 1
    return (iteration)


def Gauss_Siedel_method(A, b):
    r = [[0] for x in range(len(A))]  # r0
    res = addMatrices(multiplyMatrices(A, r), multiplyByNum(b, -1))
    iteration = 0
    norm_res = norm(res)
    while norm_res > RES_TRESHHOLD:
        if norm_res == float('inf'):
            return "no solution found"
        for i in range(len(A)):
            x = 0
            for j in range(len(A[i])):
                if j != i:
                    x += A[i][j]*r[j][0]
            r[i][0] = (b[i][0] - x)*(1/A[i][i])
        res = addMatrices(multiplyMatrices(A, r), multiplyByNum(b, -1))
        norm_res = norm(res)
        iteration += 1
    return (iteration)


def LU(A, b):
    # create L and U
    U = A
    L = IdentityMatrix(len(A))
    m = len(A)
    for k in range(m):
        for j in range(k+1, m):
            L[j][k] = U[j][k]/U[k][k]
            for i in range(k, m):
                U[j][i] -= L[j][k]*U[k][i]

    LU = addMatrices(addMatrices(L,U), multiplyByNum( IdentityMatrix(len(A)), -1))
    r = [0 for x in range(len(A))]  # r0
    for i in range(len(A)):
        sum = 0
        for j in range(i):
            sum += LU[i][j] * r[j]
        r[i] = b[i][0] - sum
    x = [0 for x in range(len(A))]
    for i in range(len(A)-1,-1,-1):
        sum = 0
        for j in range(i+1, len(A)):
            sum += LU[i][j] * x[j]
        x[i] = (1/LU[i][i])*(r[i] - sum)
    return x

def create(x):
    N = x
    b = [[sin(n * (f+1))] for n in range(N)]
    A = [[0] * N for x in range(N)]
    for i in range(N):
        A[i][i] = a1
        if i < (N-1):
            A[i][i+1] = a2
        if i >= 1:
            A[i][i - 1] = a2
        if i < (N - 2):
            A[i][i + 2] = a3
        if i >= 2:
            A[i][i - 2] = a3

    return A, b

if __name__ == '__main__':
    #initialize matrix A and vector b
    A, b = create(N)
    '''
    start = time()
    print(f'Jacobi iterations: {Jacobi_method(A, b)}')
    end = time()
    print(f'time = {end-start}')
    start = time()
    print(f'Gauss-Siedel iterations: {Gauss_Siedel_method(A, b)}')
    end = time()
    print(f'time = {end-start}')
    start = time()
    LU(A, b)
    end = time()
    print(f'LU factorization time = {end-start}')
    '''

    import matplotlib.pyplot as plt
    from matplotlib import patches

    plt.close("all")

    N_values = [100, 500, 1000, 2000, 3000]
    gauss, jacobi, lu = {0:0}, {0:0}, {0:0}
    for i in N_values:
        A, b = create(i)
        # GAUSS_SIEDEL
        start = time()
        Gauss_Siedel_method(A, b)
        end = time()
        gauss[i] = end-start
        #JACOBI
        start = time()
        Jacobi_method(A, b)
        end = time()
        jacobi[i] = end - start
        #LU
        start = time()
        #LU(A, b)
        end = time()
        lu[i] = end - start
    print(gauss.keys(), gauss.values())
    print(jacobi.keys(), jacobi.values())
    plt.plot(gauss.keys(), gauss.values(), color='g', label='Gauss-Seidel')
    plt.plot(jacobi.keys(), jacobi.values(), color='r', label='Jacobi')
    #plt.plot(lu.keys(), lu.values(), color='b', label='LU')

    plt.xlabel("N")
    plt.ylabel("Czas [s]")
    plt.title("Czas działania metod od N")
    plt.grid('both')
    plt.legend()
    plt.savefig('g-sVSj.png')
    plt.show()
