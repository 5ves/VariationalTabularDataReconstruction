from tkinter import *
from tkinter.ttk import *

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

root = Tk()

figure = Figure(figsize=(5, 4), dpi=100)
plot = figure.add_subplot(1, 1, 1)

plot.plot(0.5, 0.3, color="red", marker="o", linestyle="")

x = [ 0.1, 0.2, 0.3 ]
y = [ -0.1, -0.2, -0.3 ]
plot.plot(x, y, color="blue", marker="x", linestyle="")

canvas = FigureCanvasTkAgg(figure, root)
canvas.get_tk_widget().grid(row=0, column=0)

root.mainloop()