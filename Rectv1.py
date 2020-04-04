#
#   Версия приложения от 04.04.2020
#   Здесь не работает восстановление А. Начаты работы по выправлению формул для разрывных функций.
#   В след версии хочу убрать в интегрированиях сумммирование с пред значениями.

from mpl_toolkits import mplot3d
from scipy.optimize import *
from scipy.integrate import *
import matplotlib
import numpy as np
import random
import sys
from scipy.interpolate import griddata
import tkinter as tk
import matplotlib.pyplot as plt
#matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler

print("You are using Python {}.{}.".format(sys.version_info.major, sys.version_info.minor))

global F, ISQ, ISW, error, I, J, Tries, Z
F = 0
error = 4.0
ISW = 0
ISQ = 0
I = 0
J = 0
Tries = 100

hButton = 2
wButton = 5
hInfo = 3
wInfo = 6

m = 'SLSQP'  # 'SLSQP' #'COBYLA' # minimization method
eps = 0.01  # epsilon
lVal = 5  # значение буквы
N = 20  # grid size
IG = 0
root = tk.Tk()
NewNoise = tk.BooleanVar() # new noise control variable
NewNoise.set(1)
ForceBuild = tk.BooleanVar() #force build conttrol variable
ForceBuild.set(1)
Z = np.zeros((N - 1, N - 1))


#errorLabel = tk.Label(root, width=wInfo, height=hInfo)
errorEntry = tk.Entry(root)
moreButton = tk.Button(root, text="+ 0.1", width=wButton, height=hButton)
lessButton = tk.Button(root, text="- 0.1", width=wButton, height=hButton)
retryButton = tk.Button(root, text="Retry", width=wButton, height=hButton)
randomButton = tk.Checkbutton(text="New noise", variable=NewNoise, onvalue=1, offvalue=0)
forceBuildButton = tk.Checkbutton(text="Force build", variable=ForceBuild, onvalue=1, offvalue=0)
#errorLabel['text'] = str(error)
fig = plt.figure(dpi=70)

def Sq(u):
    global ISQ
    return (1 + u ** 2) ** (1 / 2) + ISQ


def Sw(u):
    global ISW
    return (1 + u ** 2) ** (1 / 2) + ISW


def S(u):
    global I
    global J
    global F
    row = []
    row2 = []
    row4 = []
    F[I - 1, J - 1] = ((1 + u[0] ** 2 + u[1] ** 2) ** (1 / 2))
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
            f = f + (lMatrix[i][j] * F[i][j])
    return (1 / 9) * f


def AddError():
    global error
    error = round(error + 0.1, 1)
    errorEntry.delete(0, tk.END)
    errorEntry.insert(0, error)
    #errorLabel['text'] = error
    Run()


def SubtractError():
    global error
    error = round(error - 0.1, 1)
    #errorLabel['text'] = error
    errorEntry.delete(0, tk.END)
    errorEntry.insert(0, error)
    Run()


def errorEntryValueChanged(event):
    global error
    error = round(eval(errorEntry.get()), 1)
    errorEntry.delete(0, tk.END)
    errorEntry.insert(0, error)
    Run()


def Retry():
    Run()

def FillOrigin():  # fill letter
    return FillOriginA()

def FillOriginN():  # fill letter
    Z = np.zeros((N - 1, N - 1))
    for i in range(4, N - 4):
        Z[i, 4] = lVal
        Z[i, N - 4] = lVal
        if i != 4 or i != N - 4:
            Z[i, i] = lVal
    return Z


def FillOriginA():  # fill letter
    Z = np.zeros((N - 1, N - 1))
    j = int(N/2 - 1)
    cnt = 0
    for i in range(5, N - 5):
        if cnt == 4:
            for k in range(j, j + cnt*2):
                Z[k, i] = lVal
        else:
            Z[j, i] = lVal
            Z[j + cnt * 2, i] = lVal
        cnt = cnt + 1
        j = j -1
    return Z

def AddNoise(Z):
    for z in Z:
        for i in range(len(z)):
            if z[i] == 0:
                z[i] = random.randint(2, 4)
    return Z


def RestorationRun(N, Z, e):
    global F, ISQ, ISW, error, I, J, Tries
    U = np.empty([N - 1, N - 1], float)
    Q = np.empty([N - 1, N - 1], float)
    W = np.empty([N - 1, N - 1], float)
    U[0, 0] = 0  # U00
    sq = U[0, 0]
    sw = U[0, 0]

    # filling in boundary values
    for i in range(N - 1):
        aq = Z[i, 0]
        aw = Z[0, i]
        conq = {'type': 'ineq', 'fun': lambda x: error - abs(sq + x - aq)}
        conw = {'type': 'ineq', 'fun': lambda x: error - abs(sw + x - aw)}
        resq = minimize(Sq, x0=IG, constraints=conq, method=m, options={'maxiter': Tries})
        tryCount = 0

        while resq.success == False:
            #errorEntry['bg'] = 'red'
            print("ResQ :: " + resq.message)
            #print("Resq Iterations == " + str(resq.nit) + " in " + str(i) + " " + "0")
            #print("Tries count == " + str(tryCount))
            tryCount = tryCount + 1
            if (tryCount > Tries): # and ForceBuild.get() == 0):
                return []
            resq = minimize(Sq, x0=IG, constraints=conq, method=m, options={'maxiter': Tries})

        #print("Resq Iterations == " + str(resq.nit) + " in " + str(i) + " " + str(0))
        resw = minimize(Sw, x0=IG, constraints=conw, method=m, options={'maxiter': Tries})
        tryCount = 1

        while resw.success == False:

            print("ResW (" + str(i) + " : " + str(tryCount) + ") :: " + resw.message)
            #print("Resw Iterations == " + str(resw.nit) + " in " + str(0) + " " + str(i))
            #print("Tries count == " + str(tryCount))

            tryCount = tryCount + 1
            if (tryCount > Tries):   #and ForceBuild.get() == 0):
                return []
            resw = minimize(Sw, x0=IG, constraints=conw, method=m, options={'maxiter': Tries})

        # print("Resw Iterations == " + str(resw.nit) + " in " + str(0) + " " + str(i))
        #print("conW == " + str(error - abs(sw + resw.x - aw)))
        #print("conQ == " + str(error - abs(sq + resq.x - aq)))
        if resq.x < e:
            resq.x = 0
        if resw.x < e:
            resw.x = 0

        if resq.x > 0 or resw.x > 0:
            print("jopa")

        ISW = ISW + (1 + resw.x ** 2) ** (1 / 2)
        ISQ = ISQ + (1 + resq.x ** 2) ** (1 / 2)
        Q[i, 0] = resq.x
        W[0, i] = resw.x
        U[i, 0] = resq.x + sq
        U[0, i] = resw.x + sw
        sq = sq + resq.x
        sw = sw + resw.x
        # print(res.x)
    # filling non-boundary main grid surface
    F = 0 * Z

    for i in range(1, N - 1):
        for j in range(1, N - 1):
            I = i + 1
            J = j + 1
            q = U[0, 0]
            w = U[0, 0]
            for k in range(i):
                q = q + Q[k, j]
            for k in range(j):
                w = w + W[i, k]
            a = Z[i, j]
            con = {'type': 'ineq',
                   'fun': lambda x: -abs((1 / 2) * (q + w + U[i, 0] + U[0, j] + x[0] + x[1]) - a) + error}
            res = minimize(S, x0=(IG, IG), constraints=con, method=m, options={'maxiter': Tries})
            tryCount = 1
            while res.success == False:

                print("ResMain (" + str(i) + " , " + str(j) +" : " + str(tryCount) + ") :: " + res.message)
                #print("ResMain Iterations == " + str(res.nit) + " in " + str(i) + " " + str(j))
                #print("Tries count == " + str(tryCount))

                tryCount = tryCount + 1
                if (tryCount > Tries): # and ForceBuild.get() == 0):
                    return []
                res = minimize(S, x0=(IG, IG), constraints=con, method=m, options={'maxiter': Tries})

            #print("ResMain Iterations == " + str(res.nit) + " in " + str(i) + " " + str(j))
            #print("conMain == " + str(-abs((1 / 2) * (q + w + U[i, 0] + U[0, j] + res.x[0] + res.x[1]) - a) + error))
            #print("maxcv == " + str(res.maxcv))

            if res.x[0] < eps:
                res.x[0] = 0.0
            if res.x[1] < eps:
                res.x[1] = 0.0
            Q[i, j] = res.x[0]
            W[i, j] = res.x[1]
            U[i, j] = (1 / 2) * (q + w + U[i, 0] + U[0, j] + res.x[0] + res.x[1])
            #errorEntry['bg'] = 'white'
            # q/2 + res.x[0]/2 + w/2 + res.x[1]/2 + U[i, 0] + U[0, j]
    print(error)
    return U


def Run():
    global Z
    root['bg'] = 'white'
    #errorEntry['bg'] = 'white'
    X, Y = np.meshgrid(np.linspace(0, N, N - 1), np.linspace(0, N, N - 1))  # x y grids
    O = FillOrigin()
    if NewNoise.get() == 1:
        Z = AddNoise(O.copy())
    U = RestorationRun(N, Z, eps)
    fig.clf()
    ax1 = plt.subplot(1, 3, 1, projection='3d')
    ax1.plot_surface(X, Y, O, cmap='Spectral')
    ax2 = plt.subplot(1, 3, 2, projection='3d')
    ax2.plot_surface(X, Y, Z, cmap='Spectral')
    if len(U) > 1:
        ax3 = plt.subplot(1, 3, 3, projection='3d')
        ax3.plot_surface(X, Y, U, cmap='Spectral')
    else:
        if ForceBuild.get() == 1:
            Run()
        else:
            root['bg'] = 'red'
    root.mainloop()



root.title("DATASURFACE")
root.geometry("200x150")
root.resizable(width=True, height=True)
moreButton.config(command=AddError)
lessButton.config(command=SubtractError)
retryButton.config(command=Retry)
errorEntry.bind("<Return>", errorEntryValueChanged)

errorEntry.pack()
retryButton.pack()
randomButton.pack()
forceBuildButton.pack()
moreButton.pack()
lessButton.pack()

errorEntry.delete(0, tk.END)
errorEntry.insert(0, error)

Run()
#NewNoise.set(0)
# canvas = FigureCanvasTkAgg(fig, root)
# toolbar = NavigationToolbar2Tk(canvas, root)
# toolbar.update()
#
# def on_key_press(event):
#     print("you pressed {}".format(event.key))
#     key_press_handler(event, canvas, toolbar)
#
# canvas.mpl_connect("key_press_event", on_key_press)
# canvas.get_tk_widget().pack(side="top",fill='both',expand=True)
# X = np.linspace(0, N, N-1)
# Y = np.linspace(0, N, N-1)
#ax3 = fig.add_subplot(1, 3, 1, projection='3d')
#ax1 = fig.add_subplot(1, 3, 2, projection='3d')
#ax2 = fig.add_subplot(1, 3, 3, projection='3d')
