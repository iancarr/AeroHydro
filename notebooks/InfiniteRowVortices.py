# notebook 5 - infinite row of vortices
# started 18feb2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *

# defining space
N = 100
xStart,xEnd = -6.0,6.0
yStart,yEnd = -1,1
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)              #defining the meshgrid

# defining vortex parameters
gamma = 1.0                         # defining vortex strength
nVortex = 10

# creating array of vortex position
xVortex = np.linspace(xStart,xEnd,nVortex)
yVortex = np.zeros_like(xVortex)
# xVortex,yVortex = np.meshgrid(xv,yv)

# defining a function to calculate the velocity due to a vortex
def getVelocityVortex(strength,xv,yv,X,Y):
    u = (strength/(2*pi))*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = (-strength/(2*pi))*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v
    
# defining a function to calculate vortex stream function
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = (strength/(4*pi))*(np.log((X-xv)**2+(Y-yv)**2))
    return psi

# defining the velocity matrixes to be filled
u = np.zeros_like(X)
v = np.zeros_like(Y)

# loop to produce n vortices at given positions
for i in range(nVortex):
    xv = xVortex[i]                 # defining position of current vortex
    yv = yVortex[i]
    uVortex,vVortex = getVelocityVortex(gamma,xv,yv,X,Y)
    u = u+uVortex                   # adding to final plot
    v = v+vVortex
    
# plotting finite number of vortices
size = 10
plt.figure(figsize=(size,2*(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart+2,xEnd-2)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,c='r',s=80,marker='o')
plt.show()


# INFINTE VORTEX SECTION

# defining the spacing of the vortices
a = xVortex[2]-xVortex[1]

# defining points at the center of the vortices
xInfVortex = xVortex-(a*0.5)             # not sure why these are offset by a*.5
yInfVortex = np.zeros_like(xInfVortex)

# defining shape of velocity components for infinite row
uInfVortex = np.zeros_like(u)
vInfVortex = np.zeros_like(v)

# calculating velocity field
uInfVortex = (gamma/(2*a))*(np.sinh((2*pi*Y/a))/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a)))
vInfVortex = -(gamma/(2*a))*(np.sin((2*pi*X/a))/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a)))

# plotting
size = 10
plt.figure(figsize=(size,2*(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart+2,xEnd-2)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uInfVortex,vInfVortex,density=2.0,arrowsize=1,arrowstyle='->')
plt.scatter(xInfVortex,yInfVortex,c='r',s=80,marker='o')
plt.show()