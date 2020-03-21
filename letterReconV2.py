from mpl_toolkits import mplot3d
from scipy.optimize import *
from scipy.integrate import *
import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.interpolate import griddata

def S(u):
    return (1 + u[0]**2 + u[1]**2)**(1/2)

def S1(u):
    return (1 + u**2)**(1/2)

lVal = 15

N = 20
error = 8.0

eps = 0.01
A = np.empty([N-1, N-1], float)
A.fill(np.nan)
"""
for j in range(N):
    for i in range(N):
        A[i:j] = 1/(abs(2*(i+1)/(N-1)-1) + abs(2*(j+1)/(N-1)-1)) + 1/(abs(2*(i-3)/(N-1)-1) + abs(2*(j-3)/(N-1)-1))
"""

fig = plt.figure()
X, Y = np.meshgrid(np.linspace(0, N, N-1), np.linspace(0, N, N-1))
#random.seed(version=2)
"""
Z = 1/(abs(2*(X+1)/(N-1)-1) + abs(2*(Y+1)/(N-1)-1)) + 1/(abs(2*(X-3)/(N-1)-1) + abs(2*(Y-3)/(N-1)-1)) #+ 7*(random.random() + 1)
"""
Z = 0 * X
for i in range(4, N - 4):
    Z[i, 4] = lVal
    Z[i, N-4] = lVal
    if i != 4 or i != N-4:
        Z[i, i] = lVal

ax3 = plt.subplot(1, 3, 1, projection='3d')
ax3.plot_surface(X, Y, Z, cmap='Spectral');


m = 'SLSQP'  #'SLSQP' #'COBYLA'

for z in Z:
    #random.seed(version=2)
    z += random.randint(0, 5)

ax1 = plt.subplot(1, 3, 2, projection='3d')
ax1.plot_surface(X, Y, Z, cmap='Spectral', antialiased=True);
#plt.show()

U = np.empty([N-1, N-1], float)
Q = np.empty([N-1, N-1], float)
W = np.empty([N-1, N-1], float)

U[0,0] = 0
sq = U[0,0]/2
sw = U[0,0]/2
IG = Z.max()- Z.min()
for i in range(1, N-1):
    aq = Z[i, 0]
    aw = Z[0, i]
    conq = {'type': 'ineq', 'fun': lambda x: -(abs(sq + x/2 - 2*aq) - 2*error)}
    conw = {'type': 'ineq', 'fun': lambda x: -(abs(sw + x/2 - 2*aw) - 2*error)}
    resq = minimize(S1, x0 = IG, constraints=conq, method=m)
    while resq.success == False:
        resq = minimize(S1, x0 = IG, constraints=conq, method=m)
    resw = minimize(S1, x0 = IG, constraints=conw, method=m)
    while resw.success == False:
        resw = minimize(S1, x0 = IG, constraints=conw, method=m)
    if resq.x < eps:
        resq.x = 0
    if resw.x < eps:
        resw.x = 0
    Q[i, 0] = resq.x/2
    W[0, i] = resw.x/2
    U[i, 0] = sq + resq.x
    U[0, i] = sw + resw.x
    sq = sq + resq.x/2
    sw = sw + resw.x/2
    #print(res.x)

for i in range(1, N-1):
    for j in range(1, N-1):
        q = U[0,0]/2
        w = U[0,0]/2
        for k in range(i):
            q = q + Q[k, j]/2 + Q[k, 0]/2
        for k in range(j):
            w = w + W[i, k]/2 + W[0, k]/2
        a = Z[i, j]
        con = {'type': 'ineq', 'fun': lambda x: -(abs(q + x[0]/2 + w + x[1]/2 + U[i, 0] + U[0, j] - 2*a) - 2*error)}
        res = minimize(S, x0 = (IG, IG), constraints=con, method=m)
        while res.success == False:
            res = minimize(S, x0 = (IG, IG), constraints=con, method=m)
        qlim = 0
        wlim = 0
        if res.x[0] < eps:
            res.x[0] = 0.0
        if res.x[1] < eps:
            res.x[1] = 0.0
        Q[i, j] = res.x[0]
        W[i, j] = res.x[1]
        U[i, j] = q/2 + res.x[0]/2 + w/2 + res.x[1]/2 + U[i, 0] + U[0, j]

ax2 = plt.subplot(1, 3, 3, projection='3d')
#UI = griddata((X, Y), U, (X, Y), method='cubic')

ax2.plot_surface(X, Y, U, cmap='Spectral', antialiased=True, linewidth=0)

#ax2.set_xlabel('X Label')

plt.show()
