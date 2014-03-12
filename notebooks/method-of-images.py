# notebook 7 -- Method of Images
# started 3 mar 2014
# Ian Carr

import numpy as np
import matplotlib.pyplot as plt
from math import *
from IPython.core.display import clear_output

# defining space
N = 50                      # number of points
xStart,xEnd = -2.0,2.0      # x boundaries
yStart,yEnd = -1.0,1.0      # y boundaries
x = np.linspace(xStart,xEnd,N)
y = np.linspace(xStart,xEnd,N)
X,Y = np.meshgrid(x,y)

# defining a class to compute the source attributes
class Source:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
    # get velocity field
    def velocity(self,X,Y):
        self.u = self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        self.v = self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
    # get stream function
    def streamFunction(self,X,Y):
        self.psi = self.strength/(2*pi)*np.arctan2((Y-self.y),(X-self.x))
        
# defining source parameters
strengthSource = 1.0
xSource,ySource = 0.0,0.5

# calling the class and methods within the class
source = Source(strengthSource,xSource,ySource)
source.velocity(X,Y)
source.streamFunction(X,Y)

# calling the Source class for the image of the source
sourceImage = Source(strengthSource,xSource,-ySource)
sourceImage.velocity(X,Y)
sourceImage.streamFunction(X,Y)

# superimposing the source and its image
u = source.u + sourceImage.u
v = source.v + sourceImage.v
psi = source.psi + sourceImage.psi

# plotting
size = 10
plt.figure(num=0,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=17)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(source.x,source.y,c='r',s=80,marker='o')
plt.scatter(sourceImage.x,sourceImage.y,c='r',s=80,marker='o')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4)
plt.show()

# defining class for vortex
class Vortex:
        def __init__(self,strength,x,y):
            self.strength = strength
            self.x,self.y = x,y
        # get velocity field
        def velocity(self,X,Y):
            self.u = +self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
            self.v = -self.strength/(2*pi)*(X-self.x)/((X-self.x)**2+(Y-self.y)**2)
        # get stream function
        def streamFunction(self,X,Y):
            self.psi = -self.strength/(4*pi)*np.log((X-self.x)**2+(Y-self.y)**2)
            
# defining parameters for vortex calculation
strengthVortex = 1.0
xVortex,yVortex = 0.0,0.5

# calling the class and methods for vortex
vortex = Vortex(strengthVortex,xVortex,yVortex)
vortex.velocity(X,Y)
vortex.streamFunction(X,Y)

# calling the class for the image of the vortex
vortexImage = Vortex(-strengthVortex,xVortex,-yVortex)
vortexImage.velocity(X,Y)
vortexImage.streamFunction(X,Y)

# super imposing the vortex and its image
u = vortex.u + vortexImage.u
v = vortex.v + vortexImage.v
psi = vortex.psi + vortexImage.psi

# plotting
size = 10
plt.figure(num=1,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=2.0,arrowsize=1,arrowstyle='->',linewidth=1)
plt.scatter(vortex.x,vortex.y,c='r',s=80,marker='o')
plt.scatter(vortexImage.x,vortexImage.y,c='r',s=80,marker='o')
plt.axhline(0.0,color='k',linestyle='--',linewidth=4)

# adding vortex pair to make 4 vortices total

# defining positions and strength of new vortices to be mirrored
strengthVortex = 1.0
xVortex1,yVortex1 = -0.1,0.5
xVortex2,yVortex2 = 0.1,0.5

# supplying Vortex class arg
vortex1 = Vortex(+strengthVortex,xVortex1,yVortex1)
vortex2 = Vortex(-strengthVortex,xVortex2,yVortex2)

# gradding values for Vortex methods
vortex1.velocity(X,Y)
vortex1.streamFunction(X,Y)
vortex2.velocity(X,Y)
vortex2.streamFunction(X,Y)

vortexImage1 = Vortex(-strengthVortex,xVortex1,-yVortex1)
vortexImage2 = Vortex(+strengthVortex,xVortex2,-yVortex2)

vortexImage1.velocity(X,Y)
vortexImage1.streamFunction(X,Y)
vortexImage2.velocity(X,Y)
vortexImage2.streamFunction(X,Y)

# superimposing the 4 total vortices
u = vortex1.u + vortex2.u + vortexImage1.u + vortexImage2.u
v = vortex1.v + vortex2.v + vortexImage1.v + vortexImage2.v
psi = vortex1.psi + vortex2.psi + vortexImage1.psi + vortexImage2.psi

# plotting
size = 10
plt.figure(num=2,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlim(xStart,xEnd)
plt.ylim(xStart,xEnd)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(vortex1.x,vortex1.y,s=80,c='r',marker='o')
plt.scatter(vortex2.x,vortex2.y,s=80,c='r',marker='o')
plt.scatter(vortexImage1.x,vortexImage1.y,s=80,c='r',marker='o')
plt.scatter(vortexImage2.x,vortexImage2.y,s=80,c='r',marker='o')
plt.axhline(0.0,linewidth=4,linestyle='--',c='k')

# and finally a doublet and freestream next to a wall using the method of images

Uinf = 1.0

uFreestream = Uinf*np.ones((N,N),dtype=float)
vFreestream = np.zeros_like(uFreestream)

psiFreestream = Uinf*Y

# creating class for calculation of doubelt attributes
class Doublet:
    def __init__(self,strength,x,y):
        self.strength = strength
        self.x,self.y = x,y
        
    # get doublet velocity
    def velocity(self,X,Y):
        self.u = -self.strength/(2*pi)*\
                ((X-self.x)**2-(Y-self.y)**2)/((X-self.x)**2+(Y-self.y)**2)**2
        self.v = -self.strength/(2*pi)*\
                2*(X-self.x)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)**2
                
    # get stream function
    def streamFunction(self,X,Y):
        self.psi = -self.strength/(2*pi)*(Y-self.y)/((X-self.x)**2+(Y-self.y)**2)
        
# defining doublet parameters
strengthDoublet = 1.0
xDoublet,yDoublet = 0.0,0.3

# calling Doublet class
doublet = Doublet(strengthDoublet,xDoublet,yDoublet)
doublet.velocity(X,Y)
doublet.streamFunction(X,Y)

doubletImage = Doublet(strengthDoublet,xDoublet,-yDoublet)
doubletImage.velocity(X,Y)
doubletImage.streamFunction(X,Y)

# superimposing the two vortices and freestream
u = doublet.u + doubletImage.u + uFreestream
v = doublet.v + doubletImage.v + vFreestream
psi = doublet.psi + doubletImage.psi + psiFreestream

# plotting
size = 10
plt.figure(num=3,figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(doublet.x,doublet.y,s=80,c='r',marker='o')
plt.scatter(doubletImage.x,doubletImage.y,s=80,c='r',marker='o')
plt.axhline(0.0,linewidth=4,linestyle='--',c='k')
plt.show()