import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from numpy.polynomial import Polynomial as P
x = symbols('x')
#Interpolação
##Polinomio de Lagrange
def Lagrange(pontos, valor, f):
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

def graficofLagrange(pontos, valor, f):
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
    z = np.arange(-0.9,1.5,0.01)
    
    y = []
    for i in range(len(z)):
        y.append(Pn(z[i]))

    a = []
    for i in range(len(z)):
        a.append(f.subs(x,z[i]))

    b = []
    w = []
    for i in range(len(pontos)):
        b.append(pontos[i][0])
        w.append(pontos[i][1])

    ax.plot(b,w, "r*", markersize=6, label="Pontos da tabela")
    ax.plot(z,y, label='Polinômio Interpolador P(x)')
    ax.plot(z,a, label="Função f(x)")
    ax.plot(valor,Pn(valor), "g*", markersize=6, label="Estimativa")
    ax.plot(valor,f.subs(x,valor), "yo", label="Valor exato")
    ax.legend()
    ax.grid()
    plt.show()

# def f(x): return (3+x)/(1+x)
# pontos = [[0.1, 2.82],[0.2, 2.67], [0.4, 2.43]]
# Lagrange(pontos, 1, f(x))
# graficoLagrange(pontos, 1)
# graficofLagrange(pontos, 0.25, f(x))

def Newton(pontos, valor, f):
    dif = []
    for i in range(len(pontos)):
        dif.append([])
    for i in range(len(pontos)):
        dif[0].append(pontos[i][1])
    for i in range(len(pontos)-1):
        for j in range(len(pontos)-(i+1)):
            dif[i+1].append((dif[i][j+1]-dif[i][j])/(pontos[j+i+1][0]-pontos[j][0]))
    
    Table = PrettyTable()
    points=[]
    for i in range(len(pontos)):
        points.append(pontos[i][0])
    Table.add_column("xk", points)
    for k in range(len(dif)):
        while len(dif[k]) < len(pontos):
            dif[k].append("-")
        Table.add_column("Dif_"+str(k),dif[k])

    print("Tabela")
    print(Table)

    Pn = dif[0][0]
    for i in range(1,len(dif)):
        temp = 1
        for j in range(i):
            temp *= (x-pontos[j][0])
        temp *= dif[i][0]
        Pn += temp
    
    print("Polinômio interpolador p(x) = ",end="")
    print(simplify(Pn))

    print("Polinômio interpolador avaliado em x = "+str(valor)+" é p("+str(valor)+") = ", end="")
    print(round(Pn.subs(x,valor),8))

def graficoNewton(pontos, valor):
    dif = []
    for i in range(len(pontos)):
        dif.append([])
    for i in range(len(pontos)):
        dif[0].append(pontos[i][1])
    for i in range(len(pontos)-1):
        for j in range(len(pontos)-(i+1)):
            dif[i+1].append((dif[i][j+1]-dif[i][j])/(pontos[j+i+1][0]-pontos[j][0]))
    
    Pn = dif[0][0]
    for i in range(1,len(dif)):
        temp = 1
        for j in range(i):
            temp *= (x-pontos[j][0])
        temp *= dif[i][0]
        Pn += temp
    
    fig, ax = plt.subplots()
    z = np.arange(-1,1.5,0.001)
    
    y = []
    for i in range(len(z)):
        y.append(Pn.subs(x,z[i]))

    a = []
    w = []
    for i in range(len(pontos)):
        a.append(pontos[i][0])
        w.append(pontos[i][1])

    ax.plot(z,y, label='Polinômio Interpolador P(x)')
    ax.plot(a,w, "r*", markersize=6, label="Pontos da tabela")
    ax.plot(valor,Pn.subs(x,valor), "g*", markersize=6, label="Estimativa")
    ax.legend()
    ax.grid()
    plt.show()

def graficofLagrange(pontos, valor, f):
    dif = []
    for i in range(len(pontos)):
        dif.append([])
    for i in range(len(pontos)):
        dif[0].append(pontos[i][1])
    for i in range(len(pontos)-1):
        for j in range(len(pontos)-(i+1)):
            dif[i+1].append((dif[i][j+1]-dif[i][j])/(pontos[j+i+1][0]-pontos[j][0]))
    
    Pn = dif[0][0]
    for i in range(1,len(dif)):
        temp = 1
        for j in range(i):
            temp *= (x-pontos[j][0])
        temp *= dif[i][0]
        Pn += temp
    
    fig, ax = plt.subplots()
    z = np.arange(-1,1.5,0.001)
    
    y = []
    for i in range(len(z)):
        y.append(Pn.subs(x,z[i]))

    a = []
    for i in range(len(z)):
        a.append(f.subs(x,z[i]))

    b = []
    w = []
    for i in range(len(pontos)):
        b.append(pontos[i][0])
        w.append(pontos[i][1])

    ax.plot(b,w, "r*", markersize=6, label="Pontos da tabela")
    ax.plot(z,y, label='Polinômio Interpolador P(x)')
    ax.plot(z,a, label="Função f(x)")
    ax.plot(valor,Pn.subs(x,valor), "g*", markersize=6, label="Estimativa")
    ax.plot(valor,f.subs(x,valor), "yo", label="Valor exato")
    ax.legend()
    ax.grid()
    plt.show()

# def f(x): return exp(x) + sin(x)
# pontos = [[0,1],[0.5,2.12],[1,3.55]]
# Newton(pontos,0.7,f(x))
# graficoNewton(pontos, 0.7)
# graficofLagrange(pontos, 0.7, f(x))