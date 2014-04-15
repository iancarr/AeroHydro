# lesson 11 Vortex Source Panel Method
# started 4 Mar 2013
# Ian Carr  

# importing libraries
import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# importing the airfoil geometry
coords = np.loadtxt(fname='C:/Users/Ian/Documents/GitHub/AeroHydro/resources/naca0012.dat')
xp,yp = coords[:,0],coords[:,1]

# plotting the geometry
valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,yEnd-yStart)/(xEnd-xStart)*size)
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)

#  class Panel containing the info about one panel
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                     # first point
        self.xb,self.yb = xb,yb                     # second point
        self.xc,self.yc = (xb+xa)/2,(yb+ya)/2       # center point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)   # length of the panel
        
        # orientation of a panel
        if (xb-xa<=0): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0): self.beta = pi+acos(-(yb-ya)/self.length)
        
        # location of the panel
        if (self.beta<=pi): self.loc = 'extrados'
        elif(self.beta>pi): self.loc = 'intrados'
        