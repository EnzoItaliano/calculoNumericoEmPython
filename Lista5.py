import copy
import math
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
x,t = symbols('x t')

def graficotrap(f,a,b):
    fig, ax = plt.subplots()
    z = np.arange(a,b+0.001,0.001)

    y = lambdify(x, f, "numpy")

    pontos = [[a,b],[y(a), y(b)]]

    ax.fill_between([a,b],pontos[1], color="red")
    ax.plot(z,y(z), "black")
    ax.plot(pontos[0],pontos[1])
    plt.show()


def  trapezio(f, a, b):
    h = b - a
    print("O valor de h é ", h)
    deriv = h/2*(f.subs(x,a)+f.subs(x,b))
    exact = integrate(f, (x, a, b))
    print("Pela Regra do Trapézio, temos que a integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é aproximadamente", deriv)

    print("\nComo comparação, o valor exato da integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é", exact)

    f = diff(f, x, 2)
    maior = abs(f.subs(x,a))
    if abs(f.subs(x,b)) > maior:
        maior = abs(f.subs(x,b))
    E = (-h**3/12)*maior
    print("\nLimitante")
    print("|E| <=", abs(E))

# def f(x): return log(x)+x
# trapezio(f(x), 0.5, 1)
# graficotrap(f(x), 0.5, 1)

def trapezio_gen(f, a, b, n):
    h = (b - a)/n
    print("O valor de h é", h)
    xk = np.linspace(a,b,n+1)
    fx = 0
    for i in range(len(xk)):
        if i == 0 or i == len(xk)-1:
            fx += f.subs(x, xk[i])
        else:
            fx += 2*f.subs(x, xk[i])
    deriv = (h/2)*fx
    exact = integrate(f, (x, a, b))
    print("Pela Regra do Trapézio, temos que a integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é aproximadamente", deriv)

    print("\nComo comparação, o valor exato da integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é", exact)

    f = diff(f, x, 2)
    maior = abs(f.subs(x,a))
    if abs(f.subs(x,b)) > maior:
        maior = abs(f.subs(x,b))
    E = (h**2/12)*maior*(b-a)
    print("\nLimitante")
    print("|E| <=", abs(E))


# def f(x): return sqrt(x)
# a = 1
# b = 4
# n = 6
# trapezio_gen(f(x), a, b, n)
# graficotrap(f(x),a, b)

def simpson13(f, a, b, n):
    h = (b - a)/n
    xk = np.linspace(a,b,n+1)
    fx = 0
    for i in range(len(xk)):
        if i == 0 or i == len(xk)-1:
            fx += f.subs(x, xk[i])
        elif i % 2 == 0:
            fx += 2*f.subs(x, xk[i])
        else:
            fx += 4*f.subs(x, xk[i])
    deriv = (h/3)*fx
    exact = integrate(f, (x, a, b)).evalf()

    print("Pela Regra 1/3 de Simpson, temos que a integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é aproximadamente", deriv)

    print("\nComo comparação, o valor exato da integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é", exact)

    f = diff(f, x, 4)
    maior = abs(f.subs(x,a))
    if abs(f.subs(x,b)) > maior:
        maior = abs(f.subs(x,b))
    E = (h**4/180)*maior*(b-a)
    print("\nLimitante")
    print("|E| <=", abs(E.evalf()))

def calc_parabola_vertex(x1, y1, x2, y2, x3, y3):
    denom = (x1-x2) * (x1-x3) * (x2-x3);
    A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom;
    B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom;
    C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom;
    return A,B,C

def graficoSimpson(f,a,b,n):
    xk = np.linspace(a,b,n+1)
    fx = []
    y = []
    for j in range(n-1):
        for i in range(len(xk)):
            fx.append(f.subs(x,xk[i]))
        A,B,C = calc_parabola_vertex(xk[j],fx[j],xk[j+1],fx[j+1],xk[j+2],fx[j+2])
        y.append(A*x**2 + B*x + C)

    fig, ax = plt.subplots()
    z = []
    w = []
    for i in range(n-1):
        z.append(np.arange(xk[i],xk[i+2]+0.001,0.001))

        w.append(lambdify(x, y[i], "numpy"))
        ax.plot(z[i],w[i](z[i]))

    n = np.arange(a,b+0.001,0.001)
    m = []
    for i in range(len(n)):
        m.append(f.subs(x,n[i]))
    ax.plot(n,m)
    plt.show()

def simpson_tabela13(a,b,n,y):
    h = (b - a)/n
    xk = np.linspace(a,b,n+1)
    fx = 0
    for i in range(len(xk)):
        if i == 0 or i == len(xk)-1:
            fx += y[i]
        elif i % 2 == 0:
            fx += 2*y[i]
        else:
            fx += 4*y[i]
    deriv = (h/3)*fx
    print("O valor da integral é aproximadamente", deriv)

# def f(x): return x*exp(x)+1
# a = 0
# b = 3
# n = 4
# simpson13(f(x), a, b, n)
# graficoSimpson(f(x), a, b, n)
# a = 0
# b = 6
# n = 6
# y = [0.21,0.32,0.42,0.51,0.82,0.91,1.12]
# simpson_tabela13(a,b,n,y)

def simpson38(f,a,b,n):
    h = (b-a)/n
    xk = np.linspace(a,b,n+1)
    fx = 0
    for i in range(len(xk)):
        if i == 0 or i == len(xk)-1:
            fx += f.subs(x,xk[i])
        elif i % 3 == 0:
            fx += 2*f.subs(x,xk[i])
        else:
            fx += 3*f.subs(x,xk[i])
    deriv = (3/8)*h*fx
    exact = integrate(f, (x, a, b)).evalf()

    print("Pela Regra 1/3 de Simpson, temos que a integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é aproximadamente", deriv)

    print("\nComo comparação, o valor exato da integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é", exact)

    f = diff(f, x, 4)
    maior = abs(f.subs(x,a))
    if abs(f.subs(x,b)) > maior:
        maior = abs(f.subs(x,b))
    E = (h**4/80)*maior*(b-a)
    print("\nLimitante")
    print("|E| <=", abs(E.evalf()))

def simpson_tabela38(a,b,n,y):
    h = (b - a)/n
    xk = np.linspace(a,b,n+1)
    fx = 0
    for i in range(len(xk)):
        if i == 0 or i == len(xk)-1:
            fx += y[i]
        elif i % 3 == 0:
            fx += 2*y[i]
        else:
            fx += 3*y[i]
    deriv = (3/8)*h*fx
    print("O valor da integral é aproximadamente", deriv)

# def f(x): return log(x+9)
# a = 1
# b = 7
# n = 6
# simpson38(f(x), a,b,n)
# a = 0
# b = 6
# n = 6
# y = [0.21,0.32,0.42,0.51,0.82,0.91,1.12]
# simpson_tabela38(a,b,n,y)

def quadGauss(f, a, b, n):
    table = [[[0.5773502692, 1],[-0.5773502692, 1]],
             [[0.7745966692, 0.5555555556],[0, 0.8888888889],[-0.7745966692, 0.555555556]],
             [[0.8611363116, 0.3478548451],[0.3399810436, 0.6521451549],[-0.3399810436, 0.6521451549],[0.8611363116, 0.3478548451]],
             [[0.9061798459, 0.2369268850],[0.5384693101, 0.4786286705],[0, 0.5688888889],[-0.5384693101, 0.4786286705],[-0.9061798459, 0.2369268850]]]
    
    g = f.subs(x, ((1/2)*((b-a)*t+a+b)))
    expr = lambdify(t, g, "numpy")
    if n > 2:
        soma = 0
        for i in range(n):
            soma += table[n-2][i][1]*expr(table[n-2][i][0])
        result = ((b-a)/2)*soma
    else:
        soma = 0
        for i in range(n):
            soma += expr(table[n-2][i][0])
        result = ((b-a)/2)*soma
    print("O valor aproximado da integral")
    pprint(Integral(f, (x, a, b)), use_unicode=True)
    print("é", result)

# def f(x): return x**2*log(x)
# def f(x): return exp(-x**2)
# def f(x): return sqrt(x)
# def f(x): return sin(x)/x
# def f(x): return exp(x)
# def f(x): return exp(-x**2/2)
# a = -1
# b = 1
# n = 3
# quadGauss(f(x), a, b, n)
