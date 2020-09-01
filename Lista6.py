import copy
import math
import numpy as np
from sympy import *
import plotly.graph_objects as go
from scipy.integrate import odeint
from prettytable import PrettyTable
x,y = symbols('x y')

def eulermethodf(f, a, b, x0, y0, n):
    expr = lambdify([x, y], f)
    h = (b - a) / n
    listaX = [x0]
    listaY = [y0]
    listaErr = []
    num = x0
    for i in range(n):
        num = num + h
        listaX.append(num)
    for i in range(n):
        listaY.append(listaY[len(listaY)-1] + h*expr(listaX[i],listaY[len(listaY)-1]))

    def g(y_0, xs):
        return expr(xs,y_0)

    w = odeint(g, y0,listaX)

    for i in range(len(listaY)):
        listaErr.append(abs(w[i][0] - listaY[i]))

    Table = PrettyTable()
    ks = []
    i = 0
    while(len(ks) < n + 1):
        ks.append(i)
        i+=1
    Table.add_column("i", ks)
    Table.add_column("X", listaX)
    Table.add_column("Y", listaY)
    Table.add_column("Erro", listaErr)
    print(Table)

def graficoeuler(f, a, b, x0, y0, n):
    expr = lambdify([x, y], f)
    h = (b - a) / n
    listaX = [x0]
    listaY = [y0]
    num = x0
    for i in range(n):
        num = num + h
        listaX.append(num)
    for i in range(n):
        listaY.append(listaY[len(listaY)-1] + h*expr(listaX[i],listaY[len(listaY)-1]))

    def g(y_0, xs):
        return expr(xs,y_0)

    z = np.arange(x0, listaX[len(listaX)-1]+0.001, 0.001)
    w = odeint(g, y0, z)

    c = []
    for i in range(len(w)):
        c.append(float(w[i]))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=z,y=c, name='Solução Exata'))
    fig.add_trace(go.Scatter(x=listaX, y=listaY, name="Estimativa", mode="markers"))

    fig.show()

def eulermethod(f, a, b, x0, y0, n):
    expr = lambdify([x, y], f)
    h = (b - a) / n
    listaX = [x0]
    listaY = [y0]
    num = x0
    for i in range(n):
        num = num + h
        listaX.append(num)
    for i in range(n):
        listaY.append(listaY[len(listaY)-1] + h*expr(listaX[i],listaY[len(listaY)-1]))

    Table = PrettyTable()
    ks = []
    i = 0
    while(len(ks) < n + 1):
        ks.append(i)
        i+=1
    Table.add_column("i", ks)
    Table.add_column("X", listaX)
    Table.add_column("Y", listaY)
    print(Table)

def rk2(f, a, b, x0, y0, n):
    expr = lambdify([x, y], f)
    h = (b - a) / n
    listaX = [x0]
    listaY = [y0]
    listaK1 = []
    listaK2 = []
    listaErr = []
    num = x0
    for i in range(n):
        num = num + h
        listaX.append(num)
    for i in range(n):
        listaK1.append(expr(listaX[i], listaY[i]))
        listaK2.append(expr(listaX[i] + h, listaY[i] + h*listaK1[i]))
        listaY.append(listaY[len(listaY)-1] + (h/2) * (listaK1[len(listaK1)-1] + listaK2[len(listaK2)-1]))
    listaK1.append("-")
    listaK2.append("-")

    def g(y_0, xs):
        return expr(xs,y_0)

    w = odeint(g, y0,listaX)

    for i in range(len(listaY)):
        listaErr.append(abs(w[i][0] - listaY[i]))

    Table = PrettyTable()
    ks = []
    i = 0
    while(len(ks) < n + 1):
        ks.append(i)
        i+=1
    Table.add_column("i", ks)
    Table.add_column("X", listaX)
    Table.add_column("K1", listaK1)
    Table.add_column("K2", listaK2)
    Table.add_column("Y", listaY)
    Table.add_column("Erro", listaErr)
    print(Table)

def graficork2(f, a, b, x0, y0, n):
    expr = lambdify([x, y], f)
    h = (b - a) / n
    listaX = [x0]
    listaY = [y0]
    listaK1 = []
    listaK2 = []
    num = x0
    for i in range(n):
        num = num + h
        listaX.append(num)
    for i in range(n):
        listaK1.append(expr(listaX[i], listaY[i]))
        listaK2.append(expr(listaX[i] + h, listaY[i] + h*listaK1[i]))
        listaY.append(listaY[len(listaY)-1] + (h/2) * (listaK1[len(listaK1)-1] + listaK2[len(listaK2)-1]))

    def g(y_0, xs):
        return expr(xs,y_0)

    z = np.arange(x0, listaX[len(listaX)-1]+0.001, 0.001)
    w = odeint(g, y0, z)

    c = []
    for i in range(len(w)):
        c.append(float(w[i]))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=z,y=c, name='Solução Exata'))
    fig.add_trace(go.Scatter(x=listaX, y=listaY, name="Estimativa", mode="markers"))

    fig.show()

def rk4(f, a, b, x0, y0, n):
    expr = lambdify([x, y], f)
    h = (b - a) / n
    listaX = [x0]
    listaY = [y0]
    listaK1 = []
    listaK2 = []
    listaK3 = []
    listaK4 = []
    listaErr = []
    num = x0
    for i in range(n):
        num = num + h
        listaX.append(num)
    for i in range(n):
        listaK1.append(expr(listaX[i], listaY[i]))
        listaK2.append(expr(listaX[i] + h/2, listaY[i] + (h/2)*listaK1[i]))
        listaK3.append(expr(listaX[i] + h/2, listaY[i] + (h/2)*listaK2[i]))
        listaK4.append(expr(listaX[i] + h, listaY[i] + h*listaK3[i]))
        listaY.append(listaY[len(listaY)-1] + (h/6) * (listaK1[len(listaK1)-1] + 2*listaK2[len(listaK2)-1] + 2*listaK3[len(listaK3)-1] + listaK4[len(listaK4)-1]))
    listaK1.append("-")
    listaK2.append("-")
    listaK3.append("-")
    listaK4.append("-")

    def g(y_0, xs):
        return expr(xs,y_0)

    w = odeint(g, y0,listaX)

    for i in range(len(listaY)):
        listaErr.append(abs(w[i][0] - listaY[i]))

    Table = PrettyTable()
    ks = []
    i = 0
    while(len(ks) < n + 1):
        ks.append(i)
        i+=1
    Table.add_column("i", ks)
    Table.add_column("X", listaX)
    Table.add_column("K1", listaK1)
    Table.add_column("K2", listaK2)
    Table.add_column("K3", listaK3)
    Table.add_column("K4", listaK4)
    Table.add_column("Y", listaY)
    Table.add_column("Erro", listaErr)
    print(Table)

def graficork4(f, a, b, x0, y0, n):
    expr = lambdify([x, y], f)
    h = (b - a) / n
    listaX = [x0]
    listaY = [y0]
    listaK1 = []
    listaK2 = []
    listaK3 = []
    listaK4 = []
    num = x0
    for i in range(n):
        num = num + h
        listaX.append(num)
    for i in range(n):
        listaK1.append(expr(listaX[i], listaY[i]))
        listaK2.append(expr(listaX[i] + h/2, listaY[i] + (h/2)*listaK1[i]))
        listaK3.append(expr(listaX[i] + h/2, listaY[i] + (h/2)*listaK2[i]))
        listaK4.append(expr(listaX[i] + h, listaY[i] + h*listaK3[i]))
        listaY.append(listaY[len(listaY)-1] + (h/6) * (listaK1[len(listaK1)-1] + 2*listaK2[len(listaK2)-1] + 2*listaK3[len(listaK3)-1] + listaK4[len(listaK4)-1]))
    listaK1.append("-")
    listaK2.append("-")
    listaK3.append("-")
    listaK4.append("-")

    def g(y_0, xs):
        return expr(xs, y_0)

    z = np.arange(x0, listaX[len(listaX)-1]+0.001, 0.001)
    w = odeint(g, y0, z)

    c = []
    for i in range(len(w)):
        c.append(float(w[i]))

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=z,y=c, name='Solução Exata'))
    fig.add_trace(go.Scatter(x=listaX, y=listaY, name="Estimativa", mode="markers"))

    fig.show()

# def f(x,y): return y-x
# def f(x,y): return x - y + 2
# eulermethodf(f(x,y), 0, 1, 0, 2, 4)
# graficoeuler(f(x,y), 0, 1, 0, 2, 4)
# eulermethod(f(x,y), 0, 1, 0, 2, 4)
# rk2(f(x,y), 0, 1, 0, 2, 5)
# graficork2(f(x,y), 0, 1, 0, 2, 5)
# rk4(f(x,y), 0, 1, 0, 2, 5)
# graficork4(f(x,y), 0, 1, 0, 2, 5)