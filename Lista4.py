import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial as P
x = symbols('x')
#Interpolação
##Polinomio de Lagrange
def Lagrange(pontos, valor):
    Pn = 0
    print("Polinômios coeficientes")
    for i in range(len(pontos)):
        mult = 1
        multp = 1
        div = 1
        for j in range(len(pontos)):
            if i == j: continue
            mult *= P([-pontos[j][0], 1])
            multp *= x - pontos[j][0]
            div *= pontos[i][0] - pontos[j][0]
        print("\n>>>>>>>L[%a]<<<<<<<" % i)
        pprint(multp/div)
        Pn = Pn + pontos[i][1] * (mult // div)
    print("Polinômio interpolador de Lagrange p(x) = ", end="")
    poli = list(Pn)
    for i in range(len(poli)):
        print(abs(round(poli[i],8)),end="")
        if i == 0: print(" ",end="")
        elif i == 1: print("x ", end="")
        else: print("x**%o"%i, end=" ")
        if i != len(poli)-1:
            if poli[i+1] >= 0:
                print("+ ", end="")
            else:
                print("- ", end="")
    print("\n")
    print("Polinômio interpolador avaliado em x =",valor,", é P("+str(valor)+") =" ,Pn(valor))
    
def plotL(pontos, xi, xf):
    l = []
    for i in range(len(pontos)):
        multp = 1
        div = 1
        for j in range(len(pontos)):
            if i == j: continue
            multp *= x - pontos[j][0]
            div *= pontos[i][0] - pontos[j][0]
        l.append(multp/div)

    fig, ax = plt.subplots()
    z = np.arange(xi,xf,0.01)
    y = np.zeros((len(l),len(z)))
    for i in range(len(l)):
        for j in range(len(z)):
            y[i][j] = (l[i].subs(x,z[j]))
        ax.plot(z,y[i], label=str(l[i]))
    ax.legend()
    ax.grid()
    plt.show()

def graficoLagrange(pontos, valor):
    Pn = 0
    for i in range(len(pontos)):
        mult = 1
        div = 1
        for j in range(len(pontos)):
            if i == j: continue
            mult *= P([-pontos[j][0], 1])
            div *= pontos[i][0] - pontos[j][0]
        Pn = Pn + pontos[i][1] * (mult // div)
    
    fig, ax = plt.subplots()
    z = np.arange(-4,4,0.01)
    
    y = []
    for i in range(len(z)):
        y.append(Pn(z[i]))

    x = []
    w = []
    for i in range(len(pontos)):
        x.append(pontos[i][0])
        w.append(pontos[i][1])

    ax.plot(z,y, label='Polinômio Interpolador P(x)')
    ax.plot(x,w, "r*", markersize=6, label="Pontos da tabela")
    ax.plot(valor,Pn(valor), "g*", markersize=6, label="Estimativa")
    ax.legend()
    ax.grid()
    plt.show()

pontos = [[-1, 15],[0, 8], [3, -1]]
# Lagrange(pontos, 1)
graficoLagrange(pontos, 1)