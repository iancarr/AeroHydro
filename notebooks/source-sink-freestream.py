# notebook 2 for MAE6226 started 6feb2014
# Ian Carr


import numpy as np
import matplotlib.pyplot as plt
from math import *

# defining space
N = 200                                             #num points in each direction
xStart,xEnd = -4.0,4.0                             #setting up boundaries
yStart,yEnd = -2.0,2.0
x = np.linspace(xStart,xEnd,N)                      #making the arrays
y = np.linspace(yStart,yEnd,N) 
X,Y = np.meshgrid(x,y)                              #combining the arrays into a mesh

# imposing parameters
UFs = 1.0                                            #defining freestream velocity
alphaInDegrees = 0                                  #defining angle of attack
alpha = alphaInDegrees*pi/180

# computing the velocity on grid
uFs = UFs*cos(alpha)*np.ones((N,N),dtype=float)
vFs = UFs*sin(alpha)*np.ones((N,N),dtype=float)

# computing the stream fn from velocity components
psiFs = UFs*cos(alpha)*Y - UFs*sin(alpha)*X

# defining function for velocity field of source and sink
def getVelocity(strength,xs,ys,X,Y):                #passing strength, positions, mesh
    u = (strength/(2*pi))*((X-xs)/((X-xs)**2+(Y-ys)**2))
    v = (strength/(2*pi))*((Y-ys)/((X-xs)**2+(Y-ys)**2))
    return u,v

# defining fuction to compute stream fn of source and sink
def getStreamFn(strength,xs,ys,X,Y):
    psi = (strength/(2*pi))*np.arctan((Y-ys)/(X-xs))
    return psi
    
strengthSource = 5.0                                #strength of source
xSource,ySource = -1.0,0.0

#computing the velocity components
uSource,vSource = getVelocity(strengthSource,xSource,ySource,X,Y)

#computing the stream fn
psiSource = getStreamFn(strengthSource,xSource,ySource,X,Y)

# superimposing the source and freestream
u = uFs + uSource
v = vFs + vSource
psi = psiFs + psiSource

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource,ySource,c='b',s=80,marker='o')

# computing stagnation point
xStagnation = xSource - strengthSource/(2*pi*UFs)*cos(alpha)
yStagnation = ySource - strengthSource/(2*pi*UFs)*sin(alpha)

#adding a stagnation point
plt.scatter(xStagnation,yStagnation,c='b',s=80,marker='o')

#adding a line along the stagnation streamline
if (alpha==0.0):
    plt.contour(X,Y,psi,\
        levels=[-strengthSource/2,+strengthSource/2],\
        colors='r',linewidth=2,linestyle='solid')

plt.show()
