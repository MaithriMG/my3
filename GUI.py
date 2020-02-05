import numpy as np
from scipy.io.idl import readsav
import matplotlib
from astropy.io import fits
import statistics as st
from scipy.interpolate import CubicSpline as cs
import tkinter as tk
from tkinter import *


#def show_entry_fields():
#    print("Temp: %s" %(e1.get()))#\nDwarf/Gaint: %s\nMag: %s\nseeing: %s\nExp Time: %s\nFilename: %s" % (e1.get(), e2.get(),e3.get(),e4.get(),e5.get(),e6.get()))
master = tk.Tk()
tk.Label(master,text="Temp").grid(row=0)
e1 = tk.Entry(master)
e1.grid(row=0, column=1)
tk.Button(master,text='Proceed',command=master.quit).grid(row=8,column=0, sticky=tk.W,pady=10)
master.mainloop()

RD=np.loadtxt("Combineout.Dat")
teff=RD.transpose()[0]
logg=RD.transpose()[1]
flux=RD.transpose()[2:282226]
f=(np.transpose(flux))*1e-3
RD1=np.loadtxt("3000_0.0.txt")
l=RD1.transpose()[0]
lam=l[42777:]

#Temp=int(input("Please enter the value of temperature:" ))
#G=float(input("logg value:" ))
#mag=float(input("Please enter the value of magnitude:" ))
#Seeing=float(input("At what seeing?" ))
#Exptime=float(input("The Exposure time in seconds?" ))

a=[]
for i in range(len(teff)):
    if (teff[i]==int(e1.get())):
        P=logg[i]
        a.append(P)
#print(a)
"""
a=[0,5]
root = tk.Tk()
v = tk.IntVar()
def ShowChoice():
    print(v.get())

tk.Label(root,text=Choose logg,justify = tk.LEFT,padx = 20).pack()
for val, a in enumerate(a):
    tk.Radiobutton(root,text=a,padx = 20,variable=v,command=ShowChoice,value=val).pack(anchor=tk.W)

tk.Button(root,text="next", command=root.quit).pack()
root.mainloop()
#print(v.get())
"""
#a=[0,5]
root = Tk()
frame=Frame(root,width=100, height=50)
frame.pack()

lab=Label(root,text="""Avaliable logg values""",justify = tk.LEFT,padx = 20).pack()
lab =Label(root,text=a,font=("Arial Bold", 50))
#frame.geometry('350x200')
lab.pack()
root.mainloop()

#def show():
#    print("Dwarf/Gaint(logg): %s" %(e2.get()))#\nDwarf/Gaint: %s\nMag: %s\nseeing: %s\nExp Time: %s\nFilename: %s" % (e1.get(), e2.get(),e3.get(),e4.get(),e5.get(),e6.get()))

master = tk.Tk()
tk.Label(master,text="Dwarf/Gaint(logg)").grid(row=0)
tk.Label(master,text="Magnitude").grid(row=1)
tk.Label(master,text="Seeing").grid(row=2)
tk.Label(master,text="Exp time in sec").grid(row=3)
e2 = tk.Entry(master)
e3 = tk.Entry(master)
e4 = tk.Entry(master)
e5 = tk.Entry(master)

e2.grid(row=0, column=1)
e3.grid(row=1, column=1)
e4.grid(row=2, column=1)
e5.grid(row=3, column=1)

tk.Button(master, text='proceed', command=master.quit).grid(row=8, column=0, sticky=tk.W, pady=4)
#tk.Button(master, text='Show', command=show).grid(row=8,column=1, sticky=tk.W, pady=4)
tk.mainloop()


for i in range(len(teff)):
    if (teff[i]==int(e1.get())) and (logg[i]==float(e2.get())):
        fl=f[i]
        
print(fl)

RD1=np.loadtxt("Eff_per_pixel_new.Dat")
W_Eff=RD1.transpose()[0]
E_Eff=RD1.transpose()[1]
W=W_Eff[0:133265]
E=E_Eff[0:133265]

#****************************************************************
RD=np.loadtxt("Seeing.Dat")
See=RD.transpose()[0]
per_flux_enclosed=RD.transpose()[1]
for i in range(len(See)):
    if (See[i]==float(e4.get())):
        per_flux_en=per_flux_enclosed[i]

#****************************************************************
F_new=cs(lam,fl)
F1=(F_new(W))*1e-3
star_flux=F1[66130]  #star flux emitted at 5500A
#print(W[53053],star_flux)
scaled_flux=((F1)/star_flux)*3.62e-12
scaled_flux_mag=scaled_flux*(10**(-float(e3.get())/2.5))

"""
plt.plot(W,F, "b",label="Emitted by source")
plt.legend(loc=1)
plt.xlabel("Lambda in A")
plt.ylabel("flux in J/m2/sec")
plt.show()

plt.plot(W,scaled_flux, "b",label="scaled flux")
plt.legend(loc=1)
plt.xlabel("Lambda in A")
plt.ylabel("flux in J/m2/sec")
plt.show()

plt.plot(W,scaled_flux_mag, "b",label="scaled flux")
plt.legend(loc=1)
plt.xlabel("Lambda in A")
plt.ylabel("flux in J/m2/sec")
plt.show()
"""
#*******************************************************************************************

Extinction=np.loadtxt("extinct.txt")
Ext_lambda=Extinction.transpose()[0]
Ext=Extinction.transpose()[1]

coeff10=np.polyfit(Ext_lambda,Ext,1)
f10=np.poly1d(coeff10)
Ext_lambda_new=np.linspace(Ext_lambda[0],Ext_lambda[-1],133265)
Ext_new=f10(Ext_lambda_new)
h=6.6261e-34
c=2.999e8

Ext_flux=np.zeros(len(W))
photon_num=np.zeros(len(W))
for i in range(len(W)):
    Ext_flux[i]=(10**(-(Ext_new[i]*2)/2.5))*scaled_flux_mag[i]
    photon_num[i]=(Ext_flux[i]*W[i]*1e-10*3.14*0.9/(h*c))

"""
plt.plot(W,Ext_flux, "b",label="Received at Tele")
plt.legend(loc=1)
plt.xlabel("Lambda in A")
plt.ylabel("flux in J/m2/sec")
plt.show()

plt.plot(W,photon_num, "b",label="Collected at Tele")
plt.legend(loc=1)
plt.xlabel("Lambda in A")
plt.ylabel("no.of photons")
plt.show()
"""

#**********************************************************************************************
S_N=np.zeros(len(W))
for i in range(len(E)):
    S_N[i]=(((E[i]/100)*photon_num[i]*per_flux_en*(int(e5.get())))**0.5)

"""
plt.plot(W,S_N, "b",label="S_N")
plt.legend(loc=1)
plt.xlabel("Lambda in A")
plt.ylabel("S_N for exp time=10min")
plt.show()
"""

matplotlib.use('TkAgg')
# from matplotlib import style
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


class GraphPage(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.title_label = tk.Label(self, text="SNR")
        self.title_label.pack()
        self.pack()

    def add_mpl_figure(self, fig):
        self.mpl_canvas = FigureCanvasTkAgg(fig, self)
        self.mpl_canvas.show()
        self.mpl_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2TkAgg(self.mpl_canvas, self)
        self.toolbar.update()
        self.mpl_canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)


class MPLGraph(Figure):

    def __init__(self):
        Figure.__init__(self, figsize=(5, 5), dpi=150)
        self.plot = self.add_subplot(111)
        self.plot.plot(W,S_N)


fig = MPLGraph()

root = tk.Tk()
graph_page = GraphPage(root)
graph_page.add_mpl_figure(fig)
root.mainloop()

