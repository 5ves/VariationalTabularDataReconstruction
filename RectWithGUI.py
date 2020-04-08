
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
lVal = 5  # значение буквы
global F, error, I, J, Tries, Z, cmapInd, ColorMap, U, O, Bottom, Top
Bottom = 2
Top = random.uniform(Bottom, lVal-0.01)#4.99
print("Noise lvl = " + str(Top))
F = 0
error = 4.86
I = 0
J = 0
Tries = 100

cmapArr = ['Spectral', 'Spectral_r', 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r']
cmapInd = 0

hButton = 2
wButton = 5
hInfo = 3
wInfo = 6

m = 'COBYLA'  # 'SLSQP' #'COBYLA' # minimization method
ColorMap = cmapArr[cmapInd] # 'Spectral'
eps = 0.01  # epsilon
N = 20  # grid size
IG = 0
Azim = -90
Elev = 90
Tools = tk.Tk()
recoveryRatio = tk.StringVar()
NewNoise = tk.BooleanVar() # new noise control variable
NewNoise.set(1)
ForceBuild = tk.BooleanVar() #force build conttrol variable
ForceBuild.set(1)
postProcessing = tk.BooleanVar() #force build conttrol variable
postProcessing.set(1)
randNoise = tk.BooleanVar() #force build conttrol variable
randNoise.set(1)
Z = np.zeros((N - 1, N - 1))

#errorLabel = tk.Label(root, width=wInfo, height=hInfo)
errorEntry = tk.Entry(Tools)
errlvlLabel = tk.Label(text="Error lvl")
moreButton = tk.Button(Tools, text="+ 0.1", width=wButton, height=hButton)
lessButton = tk.Button(Tools, text="- 0.1", width=wButton, height=hButton)
retryButton = tk.Button(Tools, text="Retry", width=wButton, height=hButton)
colorChangeButton = tk.Button(Tools, text="Change color", width=wButton+10, height=hButton)
colorResetButton = tk.Button(Tools, text="Reset color", width=wButton+10, height=hButton)
postProcessingButton = tk.Checkbutton(text="Post processing", variable=postProcessing, onvalue=1, offvalue=0)
randomButton = tk.Checkbutton(text="New noise", variable=NewNoise, onvalue=1, offvalue=0)
forceBuildButton = tk.Checkbutton(text="Force build", variable=ForceBuild, onvalue=1, offvalue=0)

def RandomNoiseClick():
    if randNoise.get() == 0:
        noiseEntry.config(state="normal") # readonly
        noiseEntry.config(background="white")
    else:
        noiseEntry.delete(0, tk.END)
        noiseEntry.config(background="grey")
        noiseEntry.config(state="disabled")

randNoiseButton = tk.Checkbutton(text="Randomize noise", variable=randNoise, onvalue=1, offvalue=0, command=RandomNoiseClick)
noiseEntry = tk.Entry(Tools)
noiseLabel = tk.Label(text="Noise lvl =")

recoveryRatioLabel = tk.Label(Tools, textvariable=recoveryRatio)
#errorLabel['text'] = str(error)
PlotsFig = plt.figure(dpi=70, num="PLOTS")

# filling UI
errlvlLabel.grid(row=1,column=1)
errorEntry.grid(row=1,column=2)
moreButton.grid(row=2,column=1)
lessButton.grid(row=2,column=2)
retryButton.grid(row=3, column=1)
colorChangeButton.grid(row=4,column=1)
colorResetButton.grid(row=4,column=2)
postProcessingButton.grid(row=5,column=1)
randomButton.grid(row=5,column=2)
forceBuildButton.grid(row=6,column=1)
randNoiseButton.grid(row=6,column=2)
noiseLabel.grid(row=7,column=1)
noiseEntry.grid(row=7,column=2)


def Sq(u):
    return (1 + u ** 2) ** (1 / 2)

def Sw(u):
    return (1 + u ** 2) ** (1 / 2)

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
    # for i in range(I):
        # for j in range(J):
            # f = (lMatrix[i][j] * F[i][j]) #+ f

    return (1 / 9) * (lMatrix[I-1][J-1] * F[I-1][J-1])

def AddError():
    global error
    error = round(error + 0.1, 2)
    errorEntry.delete(0, tk.END)
    errorEntry.insert(0, error)
    #errorLabel['text'] = error
    Run()

def SubtractError():
    global error
    error = round(error - 0.1, 2)
    #errorLabel['text'] = error
    errorEntry.delete(0, tk.END)
    errorEntry.insert(0, error)
    Run()

def NextColorPaint():
    global cmapInd, ColorMap
    if cmapInd == len(cmapArr) - 1:
        cmapInd = 0
        ColorMap = cmapArr[cmapInd]
    else:
        cmapInd += 1
        ColorMap = cmapArr[cmapInd]
    Run(OnlyPaint=True)

def ResetColor():
    global cmapInd, ColorMap
    cmapInd = 0
    ColorMap = cmapArr[cmapInd]
    Run(OnlyPaint=True)

def RandomNoiseClick():
    if randNoise.get() == 0:
        noiseEntry.delete(0, tk.END)
        noiseEntry.config(state=NORMAL)


def errorEntryValueChanged(event):
    global error
    error = round(eval(errorEntry.get()), 2)
    errorEntry.delete(0, tk.END)
    errorEntry.insert(0, error)
    Run()

def Retry():
    Run()

def FillOrigin():  # fill letter
    return FillOriginN() #  <-- required letter filling method here

def FillOriginN():  # N
    Z = np.zeros((N - 1, N - 1))
    for i in range(4, N - 4):
        Z[i, 4] = lVal
        Z[i, N - 4] = lVal
        if i != 4 or i != N - 4:
            Z[i, i] = lVal
    return Z

def FillOriginA():  # fill letter A
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

def FillOriginW():  # bounds test W
    Z = np.zeros((N - 1, N - 1))
    j = int(N/2 - 1)
    cnt = 0
    for i in range(0, N-1):
        Z[i, 0] = lVal
        Z[i, -1] = lVal
        if i > 9:
            Z[i, j] = lVal
            Z[i, j + 2*cnt] = lVal
            cnt = cnt + 1
            j = j - 1
    return Z

def FillOriginDK():  # D K
    Z = np.zeros((N - 1, N - 1))
    j = int(N/2 - 1)
    cnt = 0
    for i in range(4, N-4):
        Z[i, 0] = lVal
        Z[i, int(N/2 + 1)] = lVal
        if cnt == 0:
            for j in range(1, int(N/4)+1):
                Z[i, j] = lVal
            j = int(N/4)
        elif cnt <= 3:
            Z[i, j + cnt] = lVal
        else:
            Z[i, j + 3] = lVal
        Z[i, int(N/2 + 1) + 6 - cnt] = lVal
        if i > 9:
            cnt = cnt - 1
        elif i < 9:
            cnt = cnt + 1
    return Z

def AddNoise(Z):
    global Bottom, Top
    for z in Z:
        for i in range(len(z)):
            if z[i] == 0:
                z[i] += random.uniform(Bottom, Top)
    return Z

def RecoveryRatioCalculation():
    global O, U
    ratio = 0
    for i in range(N - 1):
        for j in range(N - 1):
            ratio = ratio + abs(U[i, j] - O[i, j])
    recoveryRatio.set(str(ratio))

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
            print("Resw Iterations == " + str(resw.nit) + " in " + str(0) + " " + str(i))
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

        # if resq.x > 0 or resw.x > 0:
            # print("q: " + str(resq.x) + " w: " + str(resw.x) +" i: " + str(i))

        Q[i, 0] = resq.x
        W[0, i] = resw.x

        U[i, 0] = resq.x
        U[0, i] = resw.x
        # print(res.x)
    # filling non-boundary main grid surface
    F = 0 * Z

    for i in range(1, N - 1):
        for j in range(1, N - 1):
            I = i + 1
            J = j + 1

            a = Z[i, j]
            con = {'type': 'ineq',
                   'fun': lambda x: -abs((1 / 2) * (x[0] + x[1]) - a) + error}
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

            # if res.x[0] > 0 or res.x[1] > 0:
                # print("q: " + str(res.x[0]) + " w: " + str(res.x[1]) + "  (" + str(i) + ", "+str(j)+")")

            Q[i, j] = res.x[0]
            W[i, j] = res.x[1]
            U[i, j] = (1 / 2) * (res.x[0] + res.x[1])
    print(error)
    return U

def Run(OnlyPaint = False):
    global Z, U, O
    Tools['bg'] = 'white'
    PlotsFig.clf()
    X, Y = np.meshgrid(np.linspace(0, N, N - 1), np.linspace(0, N, N - 1))  # x y grids
    if not OnlyPaint:
        O = FillOrigin()
        if NewNoise.get() == 1:
            Z = AddNoise(O.copy())
        U = RestorationRun(N, Z, eps)
        if postProcessing.get() == 1:
            for i in range(N - 1):
                for j in range(N - 1):
                    if U[i, j] > 0:
                        U[i, j] = U[i, j] + error
                        # print("(" + str(i) + ", " + str(j) + ")")
    RecoveryRatioCalculation()
    ax1 = plt.subplot(1, 3, 1, projection='3d')
    ax1.plot_surface(X, Y, O, cmap=ColorMap)
    ax2 = plt.subplot(1, 3, 2, projection='3d')
    ax2.view_init(elev = Elev, azim=Azim)
    ax1.view_init(elev = Elev, azim=Azim)
    ax2.plot_surface(X, Y, Z, cmap=ColorMap)
    if len(U) > 1:
        ax3 = plt.subplot(1, 3, 3, projection='3d')
        ax3.plot_surface(X, Y, U, cmap=ColorMap)
        ax3.view_init(elev = Elev, azim=Azim)
    else:
        if ForceBuild.get() == 1:
            Run()
        else:
            Tools['bg'] = 'red'
    Tools.mainloop()


Tools.title("DATASURFACE TOOLBOX")
Tools.geometry("300x350")
Tools.resizable(width=True, height=True)
moreButton.config(command=AddError)
lessButton.config(command=SubtractError)
retryButton.config(command=Retry)
colorChangeButton.config(command=NextColorPaint)
colorResetButton.config(command=ResetColor)
errorEntry.bind("<Return>", errorEntryValueChanged)

# errlvlLabel.pack(side="left")
# errorEntry.pack(side="right")
# retryButton.pack(side="left")
# colorChangeButton.pack()
# colorResetButton.pack()
# postProcessingButton.pack()
# randomButton.pack()
# forceBuildButton.pack()
# moreButton.pack()
# lessButton.pack()
# randNoiseButton.pack() # side="left"

# recoveryRatioLabel.pack()

errorEntry.delete(0, tk.END)
errorEntry.insert(0, error)

Run()
