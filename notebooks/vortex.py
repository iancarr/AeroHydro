# notebook 4 - vortex
# started 13feb2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *

# defining space
N = 100
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)      # defining x array
y = np.linspace(yStart,yEnd,N)      # defining y array
X,Y = np.meshgrid(x,y)              # generating grid from arrays

# defining vortex parameters
gamma = 5.0
xVortex,yVortex=0.0,0.0

# define fuction for determining vortex velocity components
def getVelocityVortex(strength,xv,yv,X,Y):
    u = (strength/(2*pi))*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = (-strength/(2*pi))*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v
    
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = (strength/(4*pi))*(np.log((X-xv)**2+(Y-yv)**2))
    return psi
    
# computing the velocity field of a vortex
uVortex,vVortex = getVelocityVortex(gamma,xVortex,yVortex,X,Y)

# computing the stream function of a vortex
psiVortex = getStreamFunctionVortex(gamma,xVortex,yVortex,X,Y)

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uVortex,vVortex,\
            density=2.0, linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,s=80,c='r',marker='o')
plt.show()

# computing the velocity of a sink
strengthSink = -1.0
xSink,ySink = 0.0,0.0

# function to compute velocity 
def getVelocitySink(strength,xs,ys,X,Y):
    u = (strength/(2*pi))*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = (strength/(2*pi))*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v
    
# function to compute stream function of a sink
def getStreamFunctionSink(strength,xs,ys,X,Y):
    psi = (strength/(2*pi))*np.arctan2((Y-ys),(X-xs))
    return psi
    
# computing the velocity of a sink
uSink,vSink = getVelocitySink(strengthSink,xSink,ySink,X,Y)

# computing the stream function of a sink
psiSink = getStreamFunctionSink(strengthSink,xSink,ySink,X,Y)

# superimposing vortex and sink
u = uVortex + uSink
v = vVortex + vSink
psi = psiVortex + psiSink

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xVortex,yVortex,c='r',s=80,marker='o')


# the following code will modify the above code to represent 
# the flow past a cylinder

# defining freestream flow
Uinf = 1.0              # velocity of the freestream flow
alpha = 0.0             # angle of attack

uFreestream = Uinf*np.cos(alpha*pi/180)*np.ones((N,N),dtype=float)
vFreestream = Uinf*np.sin(alpha*pi/180)*np.ones((N,N),dtype=float)

psiFreestream = Uinf*(np.cos(alpha*pi/180)*Y+np.sin(alpha*pi/180)*X)

# defining flow around a doublet
kappa = 1.0         # strength of the doublet
xDoublet,yDoublet = 0.0,0.0

# function to calculate flow due to doublet
def getVelocityDoublet(strength,xd,yd,X,Y):
    u = -strength/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = -strength/(2*pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    return u,v
    
# function to determine stream function of a doublet
def getStreamFunctionDoublet(strength,xd,yd,X,Y):
    psi = -strength/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi
    
# computing the velocity of a doublet
uDoublet,vDoublet = getVelocityDoublet(kappa,xDoublet,yDoublet,X,Y)

# computing the stream function of a doublet
psiDoublet = getStreamFunctionDoublet(kappa,xDoublet,yDoublet,X,Y)

# superimposing the vortex doublet and freestream flows
u = uVortex + uDoublet + uFreestream
v = vVortex + vDoublet + vFreestream
psi = psiVortex + psiDoublet + psiFreestream

# determining the stagnation points
R = sqrt(kappa/(2*pi*Uinf))

xStagn1,yStagn1 = R*cos(asin(-gamma/(4*pi*Uinf*R))),R*sin(asin(-gamma/(4*pi*Uinf*R)))
xStagn2,yStagn2 = -R*cos(asin(-gamma/(4*pi*Uinf*R))),R*sin(asin(-gamma/(4*pi*Uinf*R)))

# plotting the cylinder and stagnation points
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,\
                density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
circle = plt.Circle((0,0),radius=R,color='r',alpha=0.5)
plt.gca().add_patch(circle)
plt.scatter(xVortex,yVortex,s=80,c='r',marker='o')
plt.scatter(xStagn1,yStagn1,s=80,c='g',marker='o')
plt.scatter(xStagn2,yStagn2,s=80,c='g',marker='o')
plt.show()