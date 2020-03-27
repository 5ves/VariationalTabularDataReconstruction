from mpl_toolkits import mplot3d
from scipy.optimize import *
from scipy.integrate import *
import matplotlib.pyplot as plt
import numpy as np
import random
import sys
from scipy.interpolate import griddata
import tkinter as tk

print("You are using Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))
lVal = 15 # значение буквы

# grid size
N = 20

# error lvl
global error
error = 8.0

# epsilon
eps = 0.01

# integral functions helpers
global ISW
ISW = 0
global ISQ
ISQ = 0
global I
I = 0
global J
J = 0

root = tk.Tk()
errorLabel = tk.Label(root, width=30, height=10)
moreButton = tk.Button(root, text="+ 0.1", width=30, height=10)
lessButton = tk.Button(root, text="- 0.1", width=30, height=10)
errorLabel['text'] = str(error)

def Sq(u):
    global ISQ
    return (1 + u**2 )**(1/2) + ISQ

def Sw(u):
    global ISW
    return (1 + u**2 )**(1/2) + ISW

def S(u):
    global I
    global J
    global F
    row = []
    row2 = []
    row4 = []
    F[N-1, M-1] = ((1 + u[0]**2 + u[1]**2)**(1/2))
    for k in range(J):
        if k % 2 == 0:
            row.append(2)
            row2.append(4)
            row4.append(8)
        else:
            row.append(4)
            row2.append(8)
            row4.append(16)

    row[0] = 1
    row[-1] = 1
    row2[0] = 2
    row2[-1] = 2
    row4[0] = 4
    row4[-1] = 4

    lMatrix = [[]]

    for k in range(I):
        if k % 2 == 0:
            lMatrix.append(row4)
        else:
            lMatrix.append(row2)
    lMatrix.pop(0)
    lMatrix.insert(0, row)
    lMatrix.pop()
    lMatrix.append(row)
    f = 0
    for i in range(I):
        for j in range(J):
            f = f + (lMatrix[i, j] * F[i, j])

    return (1/9) * f

def AddError():
    global error
    error = round(error + 0.1, 1)
    errorLabel['text'] = error

def SubtractError():
    global error
    error = round(error - 0.1, 1)
    errorLabel['text'] = error


fig = plt.figure()
X, Y = np.meshgrid(np.linspace(0, N, N-1), np.linspace(0, N, N-1)) # x y grids
#X = np.linspace(0, N, N-1)
#Y = np.linspace(0, N, N-1)

# fill letter
Z = np.zeros((N-1, N-1))

for i in range(4, N - 4):
    Z[i, 4] = lVal
    Z[i, N-4] = lVal
    if i != 4 or i != N-4:
        Z[i, i] = lVal

ax3 = plt.subplot(1, 3, 1, projection='3d')
ax3.plot_surface(X, Y, Z, cmap='Spectral');

# minimization method
m = 'SLSQP'  #'SLSQP' #'COBYLA'

# add noise
for z in Z:
    z += random.randint(0, 5)

ax1 = plt.subplot(1, 3, 2, projection='3d')
ax1.plot_surface(X, Y, Z, cmap='Spectral', antialiased=True);

U = np.empty([N-1, N-1], float)
Q = np.empty([N-1, N-1], float)
W = np.empty([N-1, N-1], float)

U[0, 0] = 0  # U00
sq = U[0, 0]
sw = U[0, 0]
IG = Z.max() - Z.min()  # first step for minimization. needed max size of surface

moreButton.config(command=AddError)
lessButton.config(command=SubtractError)
errorLabel.pack()
moreButton.pack()
lessButton.pack()
#root.title("Моя первая графическая программа на Python")
#root.geometry("400x250")
#root.resizable(width=False, height=False)

#root.mainloop()


# filling in boundary values
for i in range(N-1):
    aq = Z[i, 0]
    aw = Z[0, i]
    conq = {'type': 'ineq', 'fun': lambda x: -(abs(sq + x - aq) + error)}
    conw = {'type': 'ineq', 'fun': lambda x: -(abs(sw + x - aw) + error)}
    resq = minimize(Sq, x0 = IG, constraints=conq, method=m)
    while resq.success == False:
        print("ResQ :: " + resq.message)
        resq = minimize(Sq, x0 = IG, constraints=conq, method=m)

    resw = minimize(Sw, x0 = IG, constraints=conw, method=m)
    while resw.success == False:
        print("ResW :: " + resw.message)
        resw = minimize(Sw, x0 = IG, constraints=conw, method=m)

    if resq.x < eps:
        resq.x = 0
    if resw.x < eps:
        resw.x = 0
    ISW = ISW + (1 + resw.x ** 2) ** (1 / 2)
    ISQ = ISQ + (1 + resq.x**2 )**(1/2)
    Q[i, 0] = resq.x
    W[0, i] = resw.x
    U[i, 0] = sq + resq.x
    U[0, i] = sw + resw.x
    sq = sq + resq.x
    sw = sw + resw.x
    #print(res.x)

# filling non-boundary main grid surface
global F
F = 0 * Z
for i in range(1, N-1):
    for j in range(1, N-1):
        I = i + 1
        J = j + 1
        q = U[0,0]
        w = U[0,0]
        for k in range(i):
            q = q + Q[k, j]
        for k in range(j):
            w = w + W[i, k]
        a = Z[i, j]
        #con = {'type': 'ineq', 'fun': lambda x: -(abs(q + x[0]/2 + w + x[1]/2 + U[i, 0] + U[0, j] - 2*a) - 2*error)}
        con = {'type': 'ineq',
               'fun': lambda x: -abs((1/2)*(q + w + U[i, 0] + U[0, j] + x[0] + x[1]) - a) + error}
        res = minimize(S, x0 = (IG, IG), constraints=con, method=m)
        while res.success == False:
            print(" ResMain :: " + res.message)
            res = minimize(S, x0 = (IG, IG), constraints=con, method=m)
        if res.x[0] < eps:
            res.x[0] = 0.0
        if res.x[1] < eps:
            res.x[1] = 0.0
        Q[i, j] = res.x[0]
        W[i, j] = res.x[1]
        U[i, j] = (1/2)*(q + w + U[i, 0] + U[0, j] + res.x[0] + res.x[1]) #q/2 + res.x[0]/2 + w/2 + res.x[1]/2 + U[i, 0] + U[0, j]

ax2 = plt.subplot(1, 3, 3, projection='3d')
#UI = griddata((X, Y), U, (X, Y), method='cubic')

ax2.plot_surface(X, Y, U, cmap='Spectral', antialiased=True, linewidth=0)

#ax2.set_xlabel('X Label')

plt.show()


