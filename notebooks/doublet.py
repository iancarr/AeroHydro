# notebook 3 - doublet
# started 2 feb 2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *

# defining domain
N=50
xStart,xEnd = -2.0,2.0
yStart,yEnd = -1.0,1.0
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)          # generating the mesh

# defining doublet parameters
kappa = 1.0                     # strength of the doublet
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

# plotting the streamlines associated with a doublet
size = 10
plt.figure(figsize=(size,(yStart-yEnd)/(xStart-xEnd)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uDoublet,vDoublet,\
            density=2.0,linewidth=1,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,c='r',s=80,marker='o')

# defining freestream flow
Uinf = 1.0
alpha = 0.0

# calculating the velocity of the freestream
uFreestream = Uinf*np.cos(alpha)*np.ones((N,N),dtype=float)
vFreestream = Uinf*np.sin(alpha)*np.ones((N,N),dtype=float)

# calculating the streamfunction of the doublet
psiFreestream = Uinf*(np.cos(alpha*pi/180)*Y-np.sin(alpha*pi/180)*X)

# superimposing the freestream and doublet flows
u = uFreestream + uDoublet
v = vFreestream + vDoublet
psi = psiFreestream + psiDoublet

# plotting the streamlines associated with doublet in freestream
size = 10
plt.figure(figsize=(size,(yStart-yEnd)/(xStart-xEnd)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,c='r',s=80,marker='o')
plt.contour(X,Y,psi,levels=[0.0],linewidths=2.0,colors='r',linestyles='solid')

# calculating the coefficient of pressure
Cp = 1.0-(u**2+v**2)/(Uinf**2)

# plotting Cp for doublet
size = 10
plt.figure(figsize=(size,(yStart-yEnd)/(xStart-xEnd)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
cbar.set_label('$C_p$', fontsize=16)
plt.scatter(xDoublet,yDoublet,c='r',s=80,marker='o')
plt.contour(X,Y,psi,levels=[0.0],colors='r',linewidths=2.0,linestyles='solid')
plt.show()