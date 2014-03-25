# lesson 10 - Source Panel Method
# Started 25 Mar 2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *

# data from NACA is in resources file
# reading shape data from naca file
coords = np.loadtxt(fname='../resources/naca0012.dat')
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

