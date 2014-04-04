# lesson 11 Votex Source Panel Method
# started 4 Mar 2013
# Ian Carr  

# importing libraries
import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# importing the airfoil geometry
coords = np.loadtxt(fname='C:\Users\Ian\Documents\GitHub\AeroHydro\resources')
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