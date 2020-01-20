import math
import copy
import numpy as np
import prettymatrix
from prettytable import PrettyTable
from sympy import *
# Métodos diretos
## Método de Cholesky
def Determinantes(A, ordem):
    k = 1
    tempA = [[]]
    while k <= ordem:
        for i in range(k):
            for j in range(k):
                tempA[i].append(A[i][j])
        if(round(np.linalg.det(tempA), 8) == 0): return (k,round(np.linalg.det(tempA), 8))
        if k < ordem:
            for i in range(k):
                for j in range(k):
                    tempA[i].pop()
            tempA.append([])
        k+=1
    return (0,round(np.linalg.det(tempA), 8))

def Transposta(A, ordem):
    k = ordem
    mT = []
    for c in range(k):
        mT.append([])
    for i in range(k):
        for j in range(k):
            mT[i].append( A[j][i] )
    soma = 0
    for i in range(k):
        for j in range(k):
            if(mT[i][j] != A[i][j]):
                soma += 1
    if soma == 0:
        return 1

def valoresG(A, G, ordem):
    k = ordem
    for i in range(k):
        for j in range(k):
            if(i == 0 and j == 0):
                G[0][0] = A[0][0] ** (1/2)
            elif(i == j):
                soma = 0
                for c in range(0, i):
                    soma += G[i][c]**2
                G[i][j] = (A[i][j] - soma)**(1/2)
            elif(j == 0 and i >= 1):
                G[i][0] = A[i][0]/G[0][0]
            elif(j == 1 and i >= 2):
                G[i][1] = (A[i][1] - G[i][0]*G[1][0])/G[1][1]
                G[i][1] = round(G[i][1], 8)
            elif(1 <= j and j < i):
                soma = 0
                for c in range(i):
                    soma += G[i][c]*G[j][c]
                G[i][j] = (A[i][j] - soma)/G[j][j]
    
    G = np.asarray(G)
    return G

def sistLinear(G, B, ordem):
    y = symbols('y:'+str(ordem))
    mY = []
    for i in range(len(y)):
        mY.append(y[i])
    D = np.linalg.det(G)
    tempG = G.copy()
    for j in range(ordem):
        for i in range(ordem):
            tempG[i][j] = B[i]
        tempD = np.linalg.det(tempG)
        tempG = G.copy()
        mY[j] = round(tempD/D, 8)
    mTemp = []
    for i in range(len(mY)):
        mTemp.append([mY[i]])
    mY = mTemp.copy()
    mY = np.asarray(mY)
    return mY

def Cholesky(A, B):
    ordem = np.shape(A)
    if ordem[0] == 0: return print("Isso não é uma matriz")
    if ordem[0] == ordem[1]:
        det = Determinantes(A, ordem[0])
        if( det[1] <= 0):
            print("Determinante <= 0 para k = " + str(det[0]))
            return
        trans = Transposta(A, ordem[0])
        if( not trans ):
            return print("Matriz não é igual sua transposta")

        mG = np.zeros((ordem[0],ordem[1]))
        mG = valoresG(A, mG, ordem[0])
        print(prettymatrix.matrix_to_string(mG, name='G = '))

        mY = sistLinear(mG, B, ordem[0])
        print(prettymatrix.matrix_to_string(mY, name='Y = '))

        transpG = np.zeros((ordem[0],ordem[1]))
        for i in range(ordem[0]):
            for j in range(ordem[0]):
                transpG[i][j] = mG[j][i]

        transpG = np.asarray(transpG)

        mX = sistLinear(transpG, mY, ordem[0])
        print(prettymatrix.matrix_to_string(mX, name='X = '))
        
# A = [[4,2,-4],
#      [2,10,4],
#      [-4,4,9]]
# B = [0, 
#      6, 
#      5]
# Cholesky(A,B)

## Método LU
def valoresLU(L, U, A, ordem):
    for i in range(ordem):
        for j in range(ordem):
            if i == 0:
                U[0][j] = A[0][j]
            elif i == 1 and j >= 1:
                U[1][j] =  A[1][j] - L[1][0]*U[0][j]
            elif i <= j:
                soma = 0
                for c in range(i):
                    soma += L[i][c]*U[c][j]
                U[i][j] = A[i][j] - soma
        
        for j in range(ordem):
            if i == 0 and j >= 1:
                L[j][0] = A[j][0]/U[0][0]
            elif i == 1 and j >= 2:
                L[j][1] = (A[j][1] - L[j][0]*U[0][1])/U[1][1]
            elif j > i:
                soma = 0
                for c in range(i):
                    soma += L[j][c]*U[c][i]
                L[j][i] = (A[j][i] - soma)/U[i][i]

    return L, U

def LU(A,B):
    ordem = np.shape(A)
    if ordem[0] == 0: return print("Isso não é uma matriz")
    if ordem[0] != ordem[1]: return print("Esta matriz não é quadrada")
    det = Determinantes(A, ordem[0]-1)
    if( det[1] == 0):
        return print("Determinante = 0 para k = " + str(det[0]))

    L = np.eye(ordem[0])
    U = np.zeros((ordem[0],ordem[1]))

    mLU = valoresLU(L, U, A, ordem[0])

    mY = sistLinear(L, B, ordem[0])
    mX = sistLinear(U, mY, ordem[0])

    print(prettymatrix.matrices_to_string(L, U, names=['L','U']))
    print(prettymatrix.matrix_to_string(mY, name='Y = '))
    print(prettymatrix.matrix_to_string(mX, name='X = '))

# A = [[5,2,1],[-1,4,2],[2,-3,10]]
# B = [-12, 20, 3]
# LU(A,B)

## Método de Eliminação de Gauss
def MA(A, B, ordem):
    MA = np.zeros((ordem,ordem))
    MB = np.zeros(ordem)

    for i in range(ordem):
        MB[i] = B[i]
        for j in range(ordem):
            MA[i][j] = A[i][j]

    mAum = np.zeros((ordem,ordem+1))
    for k in range(ordem-1):

        for i in range(k+1,ordem):
            MB[i] = round(B[i] - (B[k] * A[i][k])/A[k][k], 8)
            for j in range(ordem):
                MA[i][j] = round(A[i][j] - (A[k][j]*A[i][k])/A[k][k], 8)


        for i in range(ordem):
            for j in range(ordem+1):
                if j == ordem:
                    mAum[i][j] = MB[i]
                else:
                    mAum[i][j] = MA[i][j]
        mAum = np.asarray(mAum)
        print(prettymatrix.matrix_to_string(mAum, name='Matriz aumentada = '))

        for i in range(ordem):
            B[i] = MB[i]
            for j in range(ordem):
                A[i][j] = MA[i][j]
    return MA, MB

def Gauss(A, B):
    ordem = np.shape(A)
    if ordem[0] == 0: return print("Isso não é uma matriz")
    if ordem[0] != ordem[1]: return print("Esta matriz não é quadrada")
    mA, mB = MA(A, B, ordem[0])
    mX = sistLinear(mA, mB, ordem[0])
    print(prettymatrix.matrix_to_string(mX, name='X = '))

# A = [[1,2,3],
#      [3,1,0],
#      [0,3,4]]
# B = [3, 4, 3]
# Gauss(A,B)

# Métodos iterativos
## Norma Linha

def norma_linha(A):
    ordem = np.shape(A)
    if ordem[0] == 0 or ordem[1] == 0: return print("Isso não é uma matriz")
    
    somaMaior = 0
    for i in range(ordem[0]):
        soma = 0
        for j in range(ordem[0]):
            soma += abs(A[i][j])
        if soma > somaMaior:
            somaMaior = soma
        print("Soma dos elementos da linha ", i+1, " = ", round(soma, 8))
    
    return somaMaior

def norma_coluna(A):
    ordem = np.shape(A)
    if ordem[0] == 0 or ordem[1] == 0: return print("Isso não é uma matriz")
    
    somaMaior = 0
    for j in range(ordem[0]):
        soma = 0
        for i in range(ordem[0]):
            soma += abs(A[i][j])
        if soma > somaMaior:
            somaMaior = soma
        print("Soma dos elementos da coluna ", j+1, " = ", soma)
    
    return somaMaior

def norma_euclidiana(A):
    ordem = np.shape(A)
    if ordem[0] == 0 or ordem[1] == 0: return print("Isso não é uma matriz")
    soma = 0
    for i in range(ordem[0]):
        for j in range(ordem[0]):
            soma += A[i][j]**2
        
    soma = soma ** (1/2)
    return soma

def norma_inf(A):
    somaMaior = 0
    for i in range(len(A)):
        if abs(A[i]) > somaMaior:
            somaMaior = abs(A[i])
    
    return somaMaior

# A = [[3,-5,7],
#      [1,-2,4],
#      [-8,1,-7]]

# print("||A||_inf = ", norma_linha(A))
# print("||A||_1 = ", norma_coluna(A))
# print("||A||_E = ", norma_euclidiana(A))

## Método Gauss Jacobi
def condicoes(A, ordem, tipo):
    flag = 0
    for i in range(ordem):
        for j in range(ordem):
            if A[i][j] > A[i][i]:
              flag = 1
        
    F = np.zeros((ordem,ordem))
    for i in range(ordem):
        for j in range(ordem):
            if i != j:
                F[i][j] = -(A[i][j]/A[i][i])
    
    if norma_linha(F) < 1:
        return true
    elif flag == 0:
        return true
    elif tipo == 1:
        if norma_coluna(F) < 1:
            return true
    else:
        return false

def GaussJacobi(A, B, X0, e):
    ordem = np.shape(A)
    if ordem[0] == 0 or ordem[1] == 0: return print("Isso não é uma matriz")
    if condicoes(A, ordem[0], 1) == false: return print("A matriz não contém os requisitos necessários para o método")
    
    xk = []
    x_x = []
    xk.append(X0)
    x_x.append('-')
    controle = 0
    end_condition = 0
    while not end_condition:
        array = []
        for i in range(len(X0)):
            array2 = []
            temp = 0
            for j in range(len(X0)):
                if i != j:
                    temp += A[i][j]*X0[j]
            array.append(round((1/A[i][i]) * (B[i] - temp),8))

        for i in range(len(array)):
            array2.append(array[i]-xk[controle][i])
        
        xk.append(array)
        x_x.append(norma_inf(array2)/norma_inf(array))
        X0 = array
        controle += 1
        if norma_inf(array2)/norma_inf(array) < e:
            end_condition = 1
    
    Table = PrettyTable(["k", "xk", "||x^(k+1) - x^(k)||/||x^(k+1)||"])
    for k in range(0, len(xk)):
        Table.add_row([k, xk[k], x_x[k]])
    
    print(Table)

# A = [[10,2,1],[1,5,1],[2,3,10]]
# B = [14, 11, 8]
# X0 = [0,0,0]
# GaussJacobi(A,B,X0,0.01)

## Método Gauss Seidel
def Sanssenfeld(A, ordem):
    B = []
    maior = -math.inf
    for i in range(ordem):
        soma = 0
        for j in range(i):
            if i == 0: continue
            soma += abs(A[i][j]/A[i][i]) * B[j]
        temp = copy.copy(soma)
        soma = 0
        for j in range(i+1,ordem):
            soma += abs(A[i][j]/A[i][i])
        print("\u03B2", i, " = ", round(temp+soma,8))
        B.append(temp+soma)
        if B[len(B)-1] > maior:
            maior = B[len(B)-1]
    if maior < 1:
        return true
    else:
        return false


def GaussSeidel(A, B, X0, e):
    flag = 0
    ordem = np.shape(A)
    if ordem[0] == 0 or ordem[1] == 0: return print("Isso não é uma matriz")
    if condicoes(A, ordem[0], 2) == false: flag = 1
    if Sanssenfeld(A, ordem[0]) == false: flag = 2
    if flag == 1 or flag == 2: return print("A matriz não atende aos requisitos do método")

    xk = []
    x_x = []
    xk.append(X0)
    x_x.append('-')
    controle = 0
    end_condition = 0
    while not end_condition:
        array = []
        array2 = []
        for i in range(len(X0)):
            temp = 0
            for j in range(len(X0)):
                if i != j:
                    if j < i:
                        temp += A[i][j]*array[j]
                    else:
                        temp += A[i][j]*X0[j]
            array.append(round((1/A[i][i]) * (B[i] - temp),8))

        for i in range(len(array)):
            array2.append(array[i]-xk[controle][i])
        
        xk.append(array)
        x_x.append(norma_inf(array2)/norma_inf(array))
        X0 = array
        controle += 1
        if norma_inf(array2)/norma_inf(array) < e:
            end_condition = 1
    
    Table = PrettyTable(["k", "xk", "||x^(k+1) - x^(k)||/||x^(k+1)||"])
    for k in range(0, len(xk)):
        Table.add_row([k, xk[k], x_x[k]])
    
    print(Table)

# A = [[5,1,1],[3,4,1],[3,3,6]]
# B = [5, 6, 0]
# X0 = [0,0,0]
# GaussSeidel(A, B, X0, 0.01)