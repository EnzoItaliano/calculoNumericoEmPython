import copy
import numpy as np
from sympy import *
import prettymatrix
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from scipy.interpolate import interp1d
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

def NewtonGregory(pontos, valor, f):
    intervalo = pontos[1][0] - pontos[0][0]
    for i in range(1,len(pontos)):
        if pontos[i][0] - pontos[i-1][0] != intervalo:
            return print("Valores de X não são equidistantes")
    
    dif = []
    for i in range(len(pontos)):
        dif.append([])
    for i in range(len(pontos)):
        dif[0].append(pontos[i][1])
    for i in range(len(pontos)-1):
        for j in range(len(pontos)-(i+1)):
            dif[i+1].append((dif[i][j+1]-dif[i][j]))

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
        temp *= (dif[i][0]/(factorial(i)*intervalo**i))
        Pn += temp
    
    print("Polinômio interpolador p(x) = ",end="")
    print(Pn)

    print("Polinômio interpolador avaliado em x = "+str(valor)+" é p("+str(valor)+") = ", end="")
    print(round(Pn.subs(x,valor),8))    

def graficoNG(pontos, valor):
    intervalo = pontos[1][0] - pontos[0][0]
    for i in range(1,len(pontos)):
        if pontos[i][0] - pontos[i-1][0] != intervalo:
            return print("Valores de X não são equidistantes")
    
    dif = []
    for i in range(len(pontos)):
        dif.append([])
    for i in range(len(pontos)):
        dif[0].append(pontos[i][1])
    for i in range(len(pontos)-1):
        for j in range(len(pontos)-(i+1)):
            dif[i+1].append((dif[i][j+1]-dif[i][j]))
    
    Pn = dif[0][0]
    for i in range(1,len(dif)):
        temp = 1
        for j in range(i):
            temp *= (x-pontos[j][0])
        temp *= (dif[i][0]/(factorial(i)*intervalo**i))
        Pn += temp
    
    fig, ax = plt.subplots()
    z = np.arange(-1,3.5,0.001)
    
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

def graficofNG(pontos, valor, f):
    intervalo = pontos[1][0] - pontos[0][0]
    for i in range(1,len(pontos)):
        if pontos[i][0] - pontos[i-1][0] != intervalo:
            return print("Valores de X não são equidistantes")
    
    dif = []
    for i in range(len(pontos)):
        dif.append([])
    for i in range(len(pontos)):
        dif[0].append(pontos[i][1])
    for i in range(len(pontos)-1):
        for j in range(len(pontos)-(i+1)):
            dif[i+1].append((dif[i][j+1]-dif[i][j]))
    
    Pn = dif[0][0]
    for i in range(1,len(dif)):
        temp = 1
        for j in range(i):
            temp *= (x-pontos[j][0])
        temp *= (dif[i][0]/(factorial(i)*intervalo**i))
        Pn += temp
    
    fig, ax = plt.subplots()
    z = np.arange(-0.4,3,0.001)
    
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
    print(valor,Pn.subs(x,valor))
    print(valor,f.subs(x,valor))
    ax.legend()
    ax.grid()
    plt.show()

# def f(x): return x/(1+x)
# pontos = [[0,1],[1,.5],[2,2/3]]
# NewtonGregory(pontos, 1.3, f(x))
# graficoNG(pontos, 1.3)
# graficofNG(pontos, 1.3, f(x))

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

def spline(pontos, valor):
    h = []
    for i in range(1,len(pontos)):
        h.append(pontos[i][0] - pontos[i-1][0])
    # print(len(h))
    M = np.zeros((len(h)-1,len(h)-1))
    for i in range(len(h)-1):
        if i == 0:
            M[i][i]   = 2*(h[i]+h[i+1])
            M[i][i+1] = h[i+1]
        elif i == len(h)-2:
            M[i][i]   = 2*(h[i]+h[i+1])
            M[i][i-1] = h[i]
        else:
            M[i][i]   = 2*(h[i]+h[i+1])
            M[i][i-1] = h[i]
            M[i][i+1] = h[i+1]

    print(prettymatrix.matrix_to_string(M, name='Matriz = '))
    B = np.zeros((len(h)-1,1))
    for i in range(1,len(h)):
        B[i-1][0] = 6*((pontos[i+1][1]-pontos[i][1])/h[i]) - 6*((pontos[i][1]-pontos[i-1][1])/h[i-1])

    print(prettymatrix.matrix_to_string(B, name='B = '))
    mu = sistLinear(M, B, len(h)-1)
    print("Spline natural: \u03BC0 = 0, \u03BC"+str(len(h))+" = 0\n")
    print("Resolvendo o sistema linear M*Y=B, temos:")
    print('\u03BC1 = ', mu[0][0])
    print('\u03BC2 = ', mu[1][0])

    alpha = np.zeros(len(h))
    beta  = np.zeros(len(h))
    gamma = np.zeros(len(h))
    for i in range(len(h)):
        if i == 0:
            alpha[i] = ((pontos[i+1][1]-pontos[i][1])/h[i]) - ((mu[i][0]/6)*h[i]) - ((0/3)*h[i])
            beta[i] = 0/2
            gamma[i] = (mu[i][0]-0)/(6*h[i])
        elif i == len(mu):
            alpha[i] = ((pontos[i+1][1]-pontos[i][1])/h[i]) - ((0/6)*h[i]) - ((mu[i-1]/3)*h[i])
            beta[i] = mu[i-1][0]/2
            gamma[i] = (0-mu[i-1][0])/(6*h[i])
        else:
            alpha[i] = ((pontos[i+1][1]-pontos[i][1])/h[i]) - ((mu[i][0]/6)*h[i]) - ((mu[i-1]/3)*h[i])
            beta[i] = mu[i-1][0]/2
            gamma[i] = (mu[i][0]-mu[i-1][0])/(6*h[i])
    
    i = np.linspace(0,len(alpha)-1,len(alpha))
    Table = PrettyTable()
    Table.add_column("i",i)
    Table.add_column("\u03B1",alpha)
    Table.add_column("\u03B2",beta)
    Table.add_column("\u03B3",gamma)
    print("\nCoeficientes dos polinomios da spline:")
    print(Table)

    S = []
    for i in range(len(alpha)):
        # print(pontos[i][1] +(alpha[i]*(x-pontos[i][0])))
        S.append(pontos[i][1] + (alpha[i]*(x-pontos[i][0])) + (beta[i]*(x-pontos[i][0])**2) + (gamma[i]*(x-pontos[i][0])**3))
    
    print("\nSpline cúbica natural:\n")
    for i in range(len(S)):
        print("P"+str(i)+"(x) = "+str(simplify(S[i]))+" , Intervalo=["+str(pontos[i][0])+","+str(pontos[i+1][0])+"]")
    print("")

    c = 0
    for i in range(1,len(pontos)):
        intervalo = [pontos[i-1][0],pontos[i][0]]
        if valor >= intervalo[0] and valor < intervalo[1]:
            c = copy.copy(i)
            break
    print("Queremos encontrar o valor para f("+str(valor)+") então devemos usar P"+str(c-1)+" pois x = "+str(valor)+" está contido no intervalo = ",intervalo)
    print("\nLogo, a função em x = "+str(valor)+" é aproximadamente: ",S[1].subs(x,valor))


def graficoSpline(pontos, valor):
    h = []
    for i in range(1,len(pontos)):
        h.append(pontos[i][0] - pontos[i-1][0])
    # print(len(h))
    M = np.zeros((len(h)-1,len(h)-1))
    for i in range(len(h)-1):
        if i == 0:
            M[i][i]   = 2*(h[i]+h[i+1])
            M[i][i+1] = h[i+1]
        elif i == len(h)-2:
            M[i][i]   = 2*(h[i]+h[i+1])
            M[i][i-1] = h[i]
        else:
            M[i][i]   = 2*(h[i]+h[i+1])
            M[i][i-1] = h[i]
            M[i][i+1] = h[i+1]

    B = np.zeros((len(h)-1,1))
    for i in range(1,len(h)):
        B[i-1][0] = 6*((pontos[i+1][1]-pontos[i][1])/h[i]) - 6*((pontos[i][1]-pontos[i-1][1])/h[i-1])

    mu = sistLinear(M, B, len(h)-1)

    alpha = np.zeros(len(h))
    beta  = np.zeros(len(h))
    gamma = np.zeros(len(h))
    for i in range(len(h)):
        if i == 0:
            alpha[i] = ((pontos[i+1][1]-pontos[i][1])/h[i]) - ((mu[i][0]/6)*h[i]) - ((0/3)*h[i])
            beta[i] = 0/2
            gamma[i] = (mu[i][0]-0)/(6*h[i])
        elif i == len(mu):
            alpha[i] = ((pontos[i+1][1]-pontos[i][1])/h[i]) - ((0/6)*h[i]) - ((mu[i-1]/3)*h[i])
            beta[i] = mu[i-1][0]/2
            gamma[i] = (0-mu[i-1][0])/(6*h[i])
        else:
            alpha[i] = ((pontos[i+1][1]-pontos[i][1])/h[i]) - ((mu[i][0]/6)*h[i]) - ((mu[i-1]/3)*h[i])
            beta[i] = mu[i-1][0]/2
            gamma[i] = (mu[i][0]-mu[i-1][0])/(6*h[i])

    S = []
    for i in range(len(alpha)):
        # print(pontos[i][1] +(alpha[i]*(x-pontos[i][0])))
        S.append(pontos[i][1] + (alpha[i]*(x-pontos[i][0])) + (beta[i]*(x-pontos[i][0])**2) + (gamma[i]*(x-pontos[i][0])**3))

    c = 0
    for i in range(1,len(pontos)):
        intervalo = [pontos[i-1][0],pontos[i][0]]
        if valor >= intervalo[0] and valor < intervalo[1]:
            c = copy.copy(i)
            break
    
    Pn = S
    
    fig, ax = plt.subplots()

    for i in range(len(pontos)-1):
        z = np.arange(pontos[i][0],pontos[i+1][0],0.001)
        y = []
        for j in range(len(z)):
            y.append(Pn[i].subs(x,z[j]))
        ax.plot(z,y, label='Polinômio Interpolador P'+str(i)+'(x)')

    a = []
    w = []
    for i in range(len(pontos)):
        a.append(pontos[i][0])
        w.append(pontos[i][1])

    ax.plot(a,w, "r*", markersize=6, label="Pontos da tabela")
    ax.plot(valor,Pn[1].subs(x,valor), "g*", markersize=6, label="Estimativa")
    ax.legend()
    ax.grid()
    plt.show()
# pontos = [[0,3.4422],[0.5,2.2302],[1,-0.8228],[1.5,-4.6133],[2,-9.0841]]
pontos = [[3,2.5],[4.5,1],[7,2.5],[9,.5]]
# spline(pontos,5)
graficoSpline(pontos, 5)