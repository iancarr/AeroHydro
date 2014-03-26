# lesson 10 - Source Panel Method
# Started 25 Mar 2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

# data from NACA is in resources file
# reading shape data from naca file
coords = np.loadtxt(fname='/Users/ian/GitHub/AeroHydro/resources/naca0012.dat')
xp,yp = coords[0:,0],coords[:,1]            # separating data into two arrays

# plotting the geometry
# in notebook %matplotlib inline - command to show figures embedded in notebook

valX,valY = 0.1,0.2
xmin,ymax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.show()

# ------------- Airfoil panel method --------------

# defining a class to define a panel
class Panel:
    def __init__(self,xa,xb,xb,yb):
        self.xa,self.ya = xa,ya                 # first point
        self.xb,self.yb = xb,yb                 # second point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2   # control point
        self.length = sqrt((xb-xa)**2+(yb-ya)/2)# length of panel
        
        # orientation of panel
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        # location of the panel
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados'
        
        self.sigma = 0.                         # source strength
        self.vt = 0.                            # tangential velocity
        self.Cp = 0.                            # pressure coeff

# discretizing airfoil shape using circular discretization method
# defining function to perform discreization
def definePanels(N,xp,yp)
    
    R = (max(xp)-min(xp))/2                         # radius of the circle
    xc,yc = (max(xp)-min(xp))/2,(max(yp)-min(yp))/2 # center point of circle
    xCircle = xc + R*np.cos(np.linspace(0,2*pi,N+1))
    yCircle = yc + R*np.sin(np.linspace(0,2*pi,N+1))
    
    x = np.copy(xCircle[0:-1])
    y = np.empty_like(x)
    
    I = 0
    for i in range(N):
        while(I<len(xp)-1):
            if (xp[i]<=x<=[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I+= 1
        a = (yp[I+1)%len(yp)]-yp[I]/(xp[(I+1)%len(yp)]-xp[I])
        b = yp[(I+1)%len(yp)]-a*xp[(I+1)%len(xp)]
        y[i] = a*x[i]+b
        
    panel = np.empty(N,dytpe = object)
    for in in range(N):
        panel[i] = Panel(x[i],y[i],x[(i+1)%N],y[(i+1)%N])
    
    return panel
    
N = 20          # number of panels
panel = definePanels(N,xp,yp) # discretization of the geometry into panels

# plotting the geometry with the panels
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
xStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa,\
    np.append([p.ya for p in panel],panel[0].ya),\
    linestyle='-',linewidth=1,\
    marker='o',markersize=6,color='r')




