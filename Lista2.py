import math
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from sympy import *
#Raízes de Equações
##Método da Bissecção
def plot2d(f, inicio, fim):
    x = []
    while inicio <= fim+0.1:
        x.append(inicio)
        inicio += 0.1
    
    y = []
    f = f.replace("x", "x[i]")
    for i in range(len(x)):
        y.append(eval(f))
    
    plt.plot(x, y)
    plt.show()
    # print(y)

def bisseccao(f, e, a, b):
    f = f.replace("x", "y")
    
    y = a
    fa = eval(f)
    fb = eval(f)
    if fa * fb > 0:
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
            fak.append(eval(f))
            y = bk[len(bk)-1]
            fbk.append(eval(f))
            xk.append((ak[len(ak)-1] + bk[len(bk)-1])/2)
            y = xk[len(xk)-1]
            fxk.append(eval(f))
            xk_x.append('-')
        else:
            if (fak[len(fak)-1] < 0 and fxk[len(fxk)-1] < 0) or (fak[len(fak)-1] > 0 and fxk[len(fxk)-1] > 0):
                ak.append(xk[len(xk)-1])
                bk.append(bk[len(bk)-1])
            else:
                ak.append(ak[len(ak)-1])
                bk.append(xk[len(xk)-1])

            y = ak[len(ak)-1]
            fak.append(eval(f))
            y = bk[len(bk)-1]
            fbk.append(eval(f))
            xk.append((ak[len(ak)-1] + bk[len(bk)-1])/2)
            y = xk[len(xk)-1]
            fxk.append(eval(f))
            temp = xk[len(xk)-1] - xk[len(xk)-2]
            if temp < 0:
                temp = temp * -1
            xk_x.append(temp)

    Table = PrettyTable(["k", "a", "b", "f(a)", "f(b)", "x", "f(x)", "|x(k) - x(k-1)|"])
    # Table.border = False
    for k in range(times):
        Table.add_row([k, ak[k], bk[k], fak[k], fbk[k], xk[k], fxk[k], xk_x[k]])


    print(Table)

# f = "x**2-3"
# plot2d(f, 0, 2)
# bisseccao(f, 0.01, -5, -1)

## Método de Newton
def newton(f, e, a, b):
    xk = []
    xk.append(b)
    xk_x = []
    xk_x.append(0)
    end_condition = 0

    if f.subs(x,xk[len(xk)-1]) * diff(diff(f,x),x).subs(x,xk[len(xk)-1]) > 0:
        while not end_condition:
            # temp = xk[len(xk)-1] - f.subs(x,xk[len(xk)-1])
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
        # Table.border = False
        for k in range(1, len(xk)):
            Table.add_row([k, xk[k], xk_x[k]])
        
        print(Table)



x = symbols('x') #define x e y como variáveis simbólicas.
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

    Table = PrettyTable(["k", "xk", "|x(k+1) - x(k)|"])     # perguntar (caderno)
    # Table.border = False
    for k in range(2, len(xk)):
        Table.add_row([k, xk[k], xk_x[k]])
        
    print(Table)

# def f(x): return cos(x) - x
# secante(f(x), 0.00001, 0.5, 0.75)

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
    # Table.border = False
    for k in range(len(xk)):
        Table.add_row([k+2, xk[k], xk_x[k]])
        
    print(Table)

# def f(x): return cos(x)-x
# regulaFalsi(f(x), 0.00001, 0.5, math.pi/4)