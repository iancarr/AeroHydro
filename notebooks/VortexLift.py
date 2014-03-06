# notbook 6 - vortex lift
# started 25feb2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *

# defining space
N = 50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

# defining doublet parameters - code borrowed from notbook 3 doublet
kappa = 2.0                     # strength of the doublet
xDoublet,yDoublet = 0.0,0.0     # location of the doublet

# function to compute doublet velocity 
def getVelocityDoublet(strength,xd,yd,X,Y):
    u = -(strength/(2*pi))*(((X-xd)**2-(Y-yd)**2)/(((X-xd)**2+(Y-yd)**2)**2))
    v = -(strength/(2*pi))*((2*(X-xd)*(Y-yd))/(((X-xd)**2+(Y-yd)**2)**2))
    return u,v
    
# function to compute doublet stream function
def getStreamFunctionDoublet(strength,xd,yd,X,Y):
    psi = -(strength/(2*pi))*((Y-yd)/((X-xd)**2+(Y-yd)**2))
    return psi
    
# calculating velocity due to doublet
uDoublet,vDoublet = getVelocityDoublet(kappa,xDoublet,yDoublet,X,Y)

# calculating the streamfunction of a doublet
psiDoublet = getStreamFunctionDoublet(kappa,xDoublet,yDoublet,X,Y)

# defining parameters for freestream flow
Uinf = 1.0                  # freestream flow velocity

# creating velocity arrays for freestream
uFreestream = np.ones((N,N),dtype=float)
vFreestream = np.zeros_like(uFreestream)

# calculating stream function for freestream flow
psiFreestream = Uinf*Y

# superimposing the doublet and freestream flow
# which is the same as defining flow around a cylinder 
u = uDoublet + uFreestream
v = vDoublet + vFreestream
psi = psiDoublet + psiFreestream

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=2.0,arrowsize=1.0,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,s=80,c='r')

# calculating and plotting a circle over the cylinder
R = np.sqrt(kappa/(2*pi*Uinf))
circle = plt.Circle((0,0),radius=R,color='r',alpha=0.5)
plt.gca().add_patch(circle)

# calculating and adding the stagnation points
xStagn1,yStagn1 = +np.sqrt(kappa/(2*pi*Uinf)),0
xStagn2,yStagn2 = -np.sqrt(kappa/(2*pi*Uinf)),0
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],c='g',s=80,marker='o')

# defining a vortex located at the origin
gamma = 7.5            # strength of a vortex
xVortex,yVortex = 0.0,0.0

# defining functions to compute the stream fn and velcity of a vortex
# code borrowed from notebook 4 - vortex
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

# superimposing vortex, doublet, and freestream
u = uVortex + uDoublet + uFreestream
v = vVortex + vDoublet + vFreestream
psi = psiVortex + psiDoublet + psiFreestream

# plotting
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.xlabel('x',fontsize=16)
plt.xlabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=2.0,arrowsize=1.0,arrowstyle='->')

# calculating and plotting a circle over the cylinder
R = np.sqrt(kappa/(2*pi*Uinf))
circle = plt.Circle((0,0),radius=R,color='r',alpha=0.5)
plt.gca().add_patch(circle)

# calculating and plotting the stagnation points
R = np.sqrt(kappa/(2*pi*Uinf))
xStagn1,yStagn1 = +np.sqrt(R**2-(gamma/(4*pi*Uinf))**2),(-gamma/(4*pi*Uinf))
xStagn2,yStagn2 = -np.sqrt(R**2-(gamma/(4*pi*Uinf))**2),(-gamma/(4*pi*Uinf))

plt.scatter(xVortex,yVortex,s=80,c='r',marker='o')
plt.scatter([xStagn1,xStagn2],[yStagn1,yStagn2],s=80,c='g')
plt.show()

# calculating pressure coefficient on the cylinder
theta = np.linspace(0,2*pi,100)
utheta = -2*Uinf*np.sin(theta)-(gamma/(2*pi*R))
Cp = 1-(utheta**2/Uinf**2)

# pressure coefficient in the case that there is no vortex
utheta_NoVortex = -2*Uinf*np.sin(theta)
Cp_NoVortex = 1-(utheta_NoVortex**2/Uinf**2)

# plotting
size = 10
plt.figure(size,figsize=(size,size))
plt.grid(True)
plt.xlabel(r'$\theta$',fontsize=18)
plt.ylabel(r'$C_p$',fontsize=18)
plt.xlim(theta.min(),theta.max())
plt.plot(theta,Cp,color='r',linewidth=2,linestyle='-')
plt.plot(theta,Cp_NoVortex,color='g',linewidth=2,linestyle='-')
plt.legend(['vortex','no vortex'],loc='best')
plt.show()