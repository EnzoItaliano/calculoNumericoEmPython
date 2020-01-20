import math
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from sympy import *
import numpy as np
x = symbols('x')
#Raízes de Equações
##Método da Bissecção
def plot2d(f, inicio, fim):
    z = np.arange(inicio,fim,0.1)
    
    y = []
    for i in range(len(z)):
        y.append(f.subs(x,z[i]))
    
    fig, ax = plt.subplots()
    ax.set(title='Gráfico função f(x)='+str(f))
    ax.plot(z,y)
    ax.grid()
    plt.show()

def bisseccao(f, e, a, b):
    fa = f.subs(x,a)
    fb = f.subs(x,b)
    if fa * fb >= 0:
        print("Não atende ao critério f(a) * f(b) < 0")
        return
    
    k = 0
    ak = []
    bk = []
    xk = []
    fak = []
    fbk = []
    xk = []
    fxk = []
    xk_x = []
    ak.append(a)
    bk.append(b)

    kf = math.log((b-a)/e,2)-1
    times = math.ceil(kf) + 1

    for k in range(times):
        if k == 0:
            y = ak[len(ak)-1]
            fak.append(round(f.subs(x,y),9))
            y = bk[len(bk)-1]
            fbk.append(round(f.subs(x,y),9))
            xk.append((ak[len(ak)-1] + bk[len(bk)-1])/2)
            y = xk[len(xk)-1]
            fxk.append(round(f.subs(x,y),9))
            xk_x.append('-')
        else:
            if (fak[len(fak)-1] < 0 and fxk[len(fxk)-1] < 0) or (fak[len(fak)-1] > 0 and fxk[len(fxk)-1] > 0):
                ak.append(xk[len(xk)-1])
                bk.append(bk[len(bk)-1])
            else:
                ak.append(ak[len(ak)-1])
                bk.append(xk[len(xk)-1])

            y = ak[len(ak)-1]
            fak.append(round(f.subs(x,y),9))
            y = bk[len(bk)-1]
            fbk.append(round(f.subs(x,y),9))
            xk.append((ak[len(ak)-1] + bk[len(bk)-1])/2)
            y = xk[len(xk)-1]
            fxk.append(round(f.subs(x,y),9))
            temp = xk[len(xk)-1] - xk[len(xk)-2]
            if temp < 0:
                temp = temp * -1
            xk_x.append(temp)

    Table = PrettyTable(["k", "a", "b", "f(a)", "f(b)", "x", "f(x)", "|x(k) - x(k-1)|"])
    for k in range(times):
        Table.add_row([k, ak[k], bk[k], fak[k], fbk[k], xk[k], fxk[k], xk_x[k]])


    print(Table)
    print("Donde \u03B5 é aproximadamente " + str(xk[len(xk)-1]))

# def f(x): return pow(x,2)-3
# plot2d(f(x), 0, 2)
# bisseccao(f(x), 0.01, 1, 2)

## Método do Ponto Fixo
def pontoFixo(f,e,xi):
    xk = []
    xk.append(xi)
    xk_x = []
    xk_x.append("-")
    end_condition = 0
    while not end_condition:
        xk.append(f.subs(x,xk[len(xk)-1]))
        xk_x.append(abs(xk[len(xk)-1]-xk[len(xk)-2]))
        if xk_x[len(xk_x)-1] < e:
            end_condition = 1
    
    Table = PrettyTable(["k", "xk", "|x(k) - x(k-1)|"])
    for k in range(0, len(xk)):
        Table.add_row([k, xk[k], xk_x[k]])
    
    print(Table)
    print("Donde \u03B5 é aproximadamente " + str(xk[len(xk)-1]))


# def f(x): return cos(x)
# pontoFixo(f(x),10**(-2), math.pi/4)

## Método de Newton
def newton(f, e, a, b):
    xk = []
    xk.append(b)
    xk_x = []
    xk_x.append(0)
    end_condition = 0

    if f.subs(x,xk[len(xk)-1]) * diff(diff(f,x),x).subs(x,xk[len(xk)-1]) > 0:
        while not end_condition:
            func = f.subs(x,xk[len(xk)-1])
            derivate = diff(f,x).subs(x,xk[len(xk)-1])
            temp = xk[len(xk)-1] - func/derivate
            xk.append(N(temp))

            temp2 = xk[len(xk)-2] - xk[len(xk)-1]
            if temp2 < 0:
                temp2 = temp2 * -1

            xk_x.append(N(temp2))
            if xk_x[len(xk_x)-1] < e:
                end_condition = 1
            
        Table = PrettyTable(["k", "xk", "|x(k) - x(k-1)|"])
        for k in range(1, len(xk)):
            Table.add_row([k, xk[k], xk_x[k]])
        
        print(Table)

        print("Donde \u03B5 é aproximadamente " + str(xk[len(xk)-1]))

# def f(x): return x**2-2
# newton(f(x), 0.00005, 1, 2)

## Método da Secante
def secante(f, e, a, b):
    xk = []
    xk.append(a)
    xk.append(b)
    xk_x = []
    xk_x.append(0)
    xk_x.append(0)
    end_condition = 0

    while not end_condition:
        temp  = f.subs(x, xk[len(xk)-1]) * (xk[len(xk)-1] - xk[len(xk)-2])
        temp2 = f.subs(x, xk[len(xk)-1]) - f.subs(x,xk[len(xk)-2])
        temp3 = xk[len(xk)-1] - (temp/temp2)
        xk.append(temp3)

        temp4 = xk[len(xk)-1] - xk[len(xk)-2]
        
        if temp4 < 0:
            temp4 = temp4 * -1

        xk_x.append(temp4)

        if xk_x[len(xk_x)-1] < e:
            end_condition = 1

    Table = PrettyTable(["k", "xk", "|x(k+1) - x(k)|"])
    for k in range(2, len(xk)):
        Table.add_row([k, xk[k], xk_x[k]])
        
    print(Table)
    print("Donde \u03B5 é aproximadamente " + str(xk[len(xk)-1]))

# def f(x): return cos(x) - x
# secante(f(x), 10**(-5), 0.5, math.pi/4)

## Método Regula Falsi
def regulaFalsi(f, e, a, b):
    xk = []
    xk_x = []

    x0 = a
    x1 = b

    end_condition = 0

    while not end_condition:
        temp = x1 - f.subs(x, x1) * (x1 - x0) / (f.subs(x, x1) - f.subs(x, x0))

        temp2 = temp - x1
        
        if temp2 < 0:
            temp2 = temp2 * -1

        if temp2 < e:
            xk.append(temp)
            xk_x.append(temp2)
            end_condition = 1
            continue

        k = f.subs(x, temp)

        if k*f.subs(x, x1) < 0:
            x0 = x1

        x1 = temp
        xk.append(temp)
        xk_x.append(temp2)
        

    Table = PrettyTable(["k", "xk", "|x(k) - x(k-1)|"])
    for k in range(len(xk)):
        Table.add_row([k+2, xk[k], xk_x[k]])
        
    print(Table)
    print("Donde \u03B5 é aproximadamente " + str(xk[len(xk)-1]))

# def f(x): return cos(x)-x
# regulaFalsi(f(x), 0.00001, 0.5, math.pi/4)
