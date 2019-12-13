import math
from decimal import *
# Conversão de Base e Erros
## Base Binária
def dectobinDecimal(n):
    n = int(n)
    binario = ""
    while(True):
        binario = binario + str(n%2)
        n = n//2
        if n == 0:
            break
    binario = binario[::-1]
    print(binario, end="")
    return binario

# dectobinDecimal(10)

def bintodecDecimal(n):
    n = str(n)
    decimal = 0
    n = n[::-1]
    tam = len(n)
    for i in range(tam):
        if n[i] == "1":
            decimal = decimal + 2**i
    print(decimal, end="")
    return decimal

# bintodecDecimal(1010)

def bintodecFracionario(n):
    n = str(n)
    x = n.split(".")
    bintodecDecimal(x[0])
    
    temp = 0
    tam = len(x[1])
    for i in range (tam):
        if x[1][i] == "1":
            temp = temp + 2 ** -(i+1)
    temp = str(temp)
    fracionario = temp.split(".")
    print("." + fracionario[1])

# bintodecFracionario(10.1)

def dectobinFracionario(n):
    n = str(n)
    x = n.split(".")
    dectobinDecimal(x[0])

    fracionario = ""
    n = "0." + x[1]
    n = float(n)
    while(True):
        n = n * 2
        if(math.floor(n) == 1):
            fracionario = fracionario + "1"
            n = n - 1
        else:
            fracionario = fracionario + "0"
        if(n == 0):
            break
    print("." + fracionario)

# dectobinFracionario(2.5)

## Base Hexadecimal
def dectohex(n):
    n = int(n)
    hexa = ""
    while(True):
        if(n % 16 == 10):
            hexa = hexa + "A"
        elif(n % 16 == 11):
            hexa = hexa + "B"
        elif(n % 16 == 12):
            hexa = hexa + "C"
        elif(n % 16 == 13):
            hexa = hexa + "D"
        elif(n % 16 == 14):
            hexa = hexa + "E"
        elif(n % 16 == 15):
            hexa = hexa + "F"
        else:
            hexa = hexa + str(n%16)
        n = n//16
        if n == 0:
            break
    hexa = hexa[::-1]
    print(hexa, end="")
    return hexa

# dectohex(26)

def dectohexF(n):
    n = str(n)
    x = n.split(".")
    dectohex(x[0])
    hexa = ""
    n = "0." + x[1]
    n = float(n)
    while(True):
        n = n * 16
        if(math.floor(n) > 1):
            hexa = hexa + str(math.floor(n))
            n = n - math.floor(n)
        else:
            hexa = hexa + "0"
        if(n == 0):
            break
    hexa = str(hexa)
    print("." + hexa)
    return hexa

# dectohexF(10.5)

def hexstring2int(n):
    decimal = 0
    n = n[::-1]
    tam = len(n)
    for i in range(tam):
        if n[i] == "A":
            decimal = decimal + 10 * 16**i
        elif n[i] == "B":
            decimal = decimal + 11 * 16**i
        elif n[i] == "C":
            decimal = decimal + 12 * 16**i
        elif n[i] == "D":
            decimal = decimal + 13 * 16**i
        elif n[i] == "E":
            decimal = decimal + 14 * 16**i
        elif n[i] == "F":
            decimal = decimal + 15 * 16**i
        else:
            decimal = decimal + int(n[i]) * 16**i
    print(decimal)
    return decimal

# hexstring2int("FACADA")

## Aritmetica de Ponto Flutuante
def paraPontoFlut(n):
    i = 0
    while n > 1:
        n *= (10**-1)
        i+=1
    # n = str(n)
    # i = str(i)
    # print(n + " x 10^" + i)
    n = float(n)
    return n, i

def trunc(n, p):
    x, c = paraPontoFlut(n)
    c = int(c)
    trunc = ""
    x = str(x)
    for i in range(p+2): # +2 para o código desconsiderar o 0.
        trunc = trunc + x[i]
    trunc = float(trunc)
    trunc = trunc * 10**c
    w, y, z = str(trunc).partition('.')
    x = ".".join([w, z[:p]])
    print(x)

# trunc(3.1415, 3)

def arred(n, p):
    x, c = paraPontoFlut(n)

    # c = int(c)
    x = round(x, p)
    x = round(x, p)
    x = x * 10**c
    print(x)

# arred(2.667, 3)
# arred(2.664, 3)

## Erro absoluto
def erroAbs(Aex, Aaprox):
    Eabs = Aex - Aaprox
    if Eabs < 0:
        Eabs *= -1
    print(Eabs)
    return Eabs

# erroAbs(math.pi, 3.14)

## Erro relativo
def erroRel(Aex, Aaprox):
    Eabs = Aex - Aaprox
    if Eabs < 0:
        Eabs *= -1
    Erel = Eabs / Aaprox
    print(Erel)
    return Erel

# erroRel(math.pi, 3.14)

# Representação binária
def Represenbin(n):
    s = int(n[0])
    c = 0
    i, j = 1, 10
    while i < 12 and j > -1:
        x = int(n[i])
        c = c + x * (2 ** j)
        i += 1
        j -= 1

    f = 0
    i, j =  12, -1
    while i < 64 and j > -53:
        x = int(n[i])
        f = f + x * (2**j)
        i += 1
        j -= 1

    result = Decimal((-1)**s * 2**(c-1023) * (1+f))
    print(result)

Represenbin('0100000000111011100100010000000000000000000000000000000000000000')
# Represenbin('01' + 8*'0' + '11101110010001' + 40*'0')