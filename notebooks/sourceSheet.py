# notebook 8 -- Source Sheet
# started 6mar2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy import integrate

# defining space
N = 100
xStart,xEnd = -1.0,1.0
yStart,yEnd = -1.5,1.5
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)
X,Y = np.meshgrid(x,y)

# defining the velocity and stream function of a freestream
Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros_like(uFreestream)

psiFreestream = Uinf*Y

# -------------- finite source sheet -----------------

# defining a class which describes a source
class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    # get source velocity
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    # get source stream function
    def streamFunction(self,X,Y):
         self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))

# defining the sources
NSource = 11
strength = 5.0
strengthSource = strength/NSource
xSource = 0.0 
ySource = np.linspace(-1.0,1.0,NSource)

# creating source line
source = np.empty(NSource,dtype=object)
for i in range(NSource):
    source[i] = Source(strengthSource,xSource,ySource[i])
    source[i].velocity(X,Y)

# superposition 
u = uFreestream.copy()
v = vFreestream.copy()

for s in source:
    u = np.add(u,s.u)
    v = np.add(v,s.v)
    
# plotting
size = 8
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xSource*np.ones(NSource,dtype=float),ySource,c='r',s=80,marker='o')
cont = plt.contourf(X,Y,np.sqrt(u**2+v**2),levels=np.linspace(0.0,0.1,10))
cbar =  plt.colorbar(cont)
cbar.set_label('U',fontsize=16);

# -------------- infinite source sheet -----------------

# the following function has been imported in scipy library
# defining the limits of the integral for velocity and stream function

a = 1.0             # one end of the integral
f = lambda x : x+a
print 'x+a =',f(1)

print integrate.quad(lambda x:x, 0.0,1.0)[0]

# defining source sheet parameters
sigma = 2.5 # strength of the source sheet

uPanel = np.empty((N,N),dtype=float)
vPanel = np.empty((N,N),dtype=float)

# boundaries of the sheet
ymin,ymax = -1.0,1.0

# computing the  velocity field
for i in range(N):
    for j in range(N):
        func = lambda s : X[i,j]/(X[i,j]**2+(Y[i,j]-s)**2)
        uPanel[i,j] = sigma/(2*pi)*integrate.quad(func,ymin,ymax)[0]
        
        func = lambda s : (Y[i,j]-s)/(X[i,j]**2+(Y[i,j]-s)**2)
        vPanel[i,j] = sigma/(2*pi)*integrate.quad(func,ymin,ymax)[0]
        
# superposition with the uniform flow field
u2 = np.add(uFreestream,uPanel)
v2 = np.add(vFreestream,vPanel)

# plotting
size = 8
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u2,v2,density=2,arrowsize=1,arrowstyle='->')
plt.axvline(0.0,(ymin-yStart)/(yEnd-yStart),(ymax-yStart)/(yEnd-yStart),\
            color='r',linewidth=3)
cont = plt.contourf(X,Y,np.sqrt(u2**2+v2**2),levels=np.linspace(0.0,0.1,10))
cbar = plt.colorbar(cont)
cbar.set_label('U',fontsize=16)
plt.show()