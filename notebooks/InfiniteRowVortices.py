# notebook 5 - infinite row of vortices
# started 18feb2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *

# defining space
N = 100
xStart,xEnd = -2.0,2.0
yStart,yEnd = -0.5,0.5
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)              #defining the meshgrid

# defining vortex parameters
gamma = 1.0                         # defining vortex strength
nVortex = 10

# creating array of vortex position
xv = np.linspace(xStart,xEnd,nVortex)
yv = np.zeros_like(xv)
xVortex,yVortex = np.meshgrid(xv,yv)

# defining a function to calculate the velocity due to a vortex
def getVortexVelocity(strength,xv,yv,X,Y):
    u = (strength/(2*pi))*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = (-strength/(2*pi))*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v
    
# defining a function to calculate vortex stream function
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = (strength/(4*pi))*(np.log((X-xv)**2+(Y-yv)**2))
    return psi
    
# looping through position and defining sheet of vortices

groups = [self.getGroup(i,header+i) for i in range(3)]
    
    