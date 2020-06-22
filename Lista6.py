import copy
import math
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
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


def f(x,y): return y-x
eulermethodf(f(x,y), 0, 1, 0, 2, 4)